#include "FoamYadeMPI.H"
#include "PstreamGlobals.H" 
#include <mpi.h> 
#include <cmath> 
#include <string>

void Foam::FoamYadeMPI::setScalarProperties(scalar _rhoP, scalar _rhoF, scalar _nu){
	rhoP = _rhoP; rhoF = _rhoF; nu = _nu; 
}

void Foam::FoamYadeMPI::printMsg(const std::string& msg){
	std::cout << " Rank : " << worldRank << "  " << msg << std::endl;   
}


void Foam::FoamYadeMPI::getRankSize(){
	if (!rankSizeSet){
		// get local rank and size 
		MPI_Comm_rank(PstreamGlobals::MPI_COMM_FOAM, &localRank); 
		MPI_Comm_size(PstreamGlobals::MPI_COMM_FOAM, &localCommSize);
		// world comm and size 
		MPI_Comm_rank(MPI_COMM_WORLD, &worldRank); 
		MPI_Comm_size(MPI_COMM_WORLD, &worldCommSize); 
		
		// diff in comm size 
		commSzDff = abs(worldCommSize-localCommSize); 
		rankSizeSet = true; 
		// alloc vector of yadeProcs, yade Master does not have participate in communications. 
		yadeProcs.resize(commSzDff-1);
		// we do not receive from yade Master, 
		for (unsigned int i=0; i != yadeProcs.size(); ++i) {
			yadeProcs[i].yRank = (int)i+1; 
			yadeProcs[i].numParticlesProc.reserve(localCommSize); 
		}
		mshTree.build_tree();  // build tree. 
		sendMeshBbox();
	}
}


void Foam::FoamYadeMPI::initFields(){
	alpha = 1.0; 
	forAll(uSource,cellI){
		
		uSource[cellI].x() = small; 
		uSource[cellI].y() = small; 
		uSource[cellI].z() = small; 
		
		uParticle[cellI].x() = small;
		uParticle[cellI].y() = small; 
		uParticle[cellI].z() = small;
		
		uSourceDrag[cellI] = small; 
	}
	//init sigmaInterp, sigmaPi & interpRangeCu
	interpRange = 2*std::pow(mesh.V()[0], 1.0/3.0);
	sigmaInterp = interpRange*0.42460; // interp_range/(2sqrt(2ln(2))) filter width half maximum;
	interpRangeCu = std::pow(interpRange, 3.0); 
	sigmaPi = 1.0/(std::pow(2*M_PI*sigmaInterp*sigmaInterp, 1.5));
}



void Foam::FoamYadeMPI::sendMeshBbox(){
	//send mesh bbox from point field, min 
	if (!rankSizeSet) {std::cerr << " rank and commsize not set. " << std::endl; return; }
	// get the mesh pointfield (cell vertices), 
	
	point minBound(1e+50, 1e+50, 1e+50); 
	point maxBound(-1e+50, -1e+50, -1e+50);  
	
	
	for (const auto& pt : mesh.points()){
		minBound.x() = Foam::min(pt.x(), minBound.x());
		minBound.y() = Foam::min(pt.y(), minBound.y());
		minBound.z() = Foam::min(pt.z(), minBound.z());
		
		maxBound.x() = Foam::max(pt.x(), maxBound.x());
		maxBound.y() = Foam::max(pt.y(), maxBound.y());
		maxBound.z() = Foam::max(pt.z(), maxBound.z());
	}
	
	std::vector<double> meshBbox = {minBound.x(), minBound.y(), minBound.z(), maxBound.x(), maxBound.y(), maxBound.z()}; 
	
	// send bounding box to every yade proc including yade Master. 
	for  (int rnk=0; rnk != commSzDff; ++rnk){
		MPI_Request req; 
		MPI_Isend(&meshBbox.front(), 6, MPI_DOUBLE, rnk, TAG_GRID_BBOX, MPI_COMM_WORLD, &req); 
		reqVec.push_back(req);
	}
	
	for (auto rq : reqVec){
		MPI_Status status; 
		MPI_Wait(&rq, &status);
	}
	
	reqVec.clear(); 
	initFields();
}


void Foam::FoamYadeMPI::recvYadeIntrs(){
  
	for (int rnk = 0; rnk != commSzDff-1; ++rnk){
		MPI_Status status; 
		MPI_Recv(&yadeProcs[rnk].numParticlesProc.front(), localCommSize, MPI_INT, rnk+1, TAG_SZ_BUFF, MPI_COMM_WORLD, &status); 
	}
	
	// those yade procs intersecting current grid. 
	for (auto& yProc : yadeProcs) {
		if (yProc.numParticlesProc[localRank] > 0 ){
			yProc.numParticles = yProc.numParticlesProc[localRank]; 
			if (!fibreCpl) {
				yProc.particleDataBuff.resize(yProc.numParticles*10);
			} 
			else {
				yProc.particleDataBuff.resize(yProc.numParticles*15);
			}
			yProc.foundBuff.resize(yProc.numParticles); // allocate 'foundbuff' for mesh search result, 
			yProc.hydroForceBuff.resize(6*yProc.numParticles); 
			std::fill(yProc.hydroForceBuff.begin(), yProc.hydroForceBuff.end(), 0.0); // init zero force here. 
			inCommProcs.push_back(std::make_shared<YadeProc>(yProc));  // keep list of yadeProcs with intersection. 
		}
	}
	
	
	// get the data from intersecting procs. 
	for (auto& yProc : inCommProcs) {
		MPI_Status status; 
		std::vector<double>& dBuff = yProc->particleDataBuff; 
		MPI_Recv(&dBuff.front(), (int) dBuff.size(), MPI_DOUBLE, yProc->yRank, TAG_YADE_DATA, MPI_COMM_WORLD, &status); 
	}

}


void Foam::FoamYadeMPI::locateAllParticles(){
	for (auto& yProc : inCommProcs){
		for (int np = 0; np != yProc->numParticles; ++np) {
			std::shared_ptr<YadeParticle> yParticle = std::make_shared<YadeParticle>(); 
			if (!fibreCpl){
				yParticle->pos.x() = yProc->particleDataBuff[np*10]; 
				yParticle->pos.y() = yProc->particleDataBuff[np*10+1]; 
				yParticle->pos.z() = yProc->particleDataBuff[np*10+2]; 
			} else {
				yParticle->pos.x() = yProc->particleDataBuff[np*15];
				yParticle->pos.y() = yProc->particleDataBuff[np*15+1];
				yParticle->pos.z() = yProc->particleDataBuff[np*15+2]; 
			}
			if (locatePt(yParticle)){
				
				yParticle->indx = np; // keep track of 'found particle' index 
				yParticle->inProc = worldRank; // redundant
				yParticle->linearVelocity.x() = yProc->particleDataBuff[np*10+3]; 
				yParticle->linearVelocity.y() = yProc->particleDataBuff[np*10+4]; 
				yParticle->linearVelocity.z() = yProc->particleDataBuff[np*10+5]; 
				
				yParticle->rotationalVelocity.x() = yProc->particleDataBuff[np*10+6]; 
				yParticle->rotationalVelocity.y() = yProc->particleDataBuff[np*10+7]; 
				yParticle->rotationalVelocity.z() = yProc->particleDataBuff[np*10+8]; 
				
				yParticle->dia = 2*yProc->particleDataBuff[np*10+9]; 
				yParticle->calcPartVol(); 
				yProc->foundBuff[np] = 1; 
				yProc->foundParticles.push_back(yParticle); 
			} 
		} 
	}
	
	// put this in a separate function? 
	
	// send the 'foundBuff' of intersecting procs to Yade. 
	for (const auto& yProc : inCommProcs) {
		std::vector<int>& fBuff = yProc->foundBuff; 
		int sz =  int (fBuff.size()); 
		MPI_Send(&fBuff.front(), sz, MPI_INT, yProc->yRank, TAG_SEARCH_RES, MPI_COMM_WORLD);  
	}
	
}


bool Foam::FoamYadeMPI::locatePt(const std::shared_ptr<YadeParticle>& yParticle){
	//TODO :: find thos particles which are in processorBoundary faces/cells. 
	bool value = false;
	
	//std::string st(std::to_string(yParticle->pos.x()) + " " + std::to_string(yParticle->pos.y())+ " " + std::to_string(yParticle->pos.z())); 
	//printMsg("searching  for " + st);
	if (isGaussianInterp) {
		yParticle->cellIds = mshTree.nnearestCellsRange(yParticle->pos, interpRange, isGaussianInterp); 
		if (yParticle->cellIds.size()) value = true; 
	}else {
		yParticle->inCell = mesh.findCell(yParticle->pos); 
		if (yParticle->inCell > -1) value = true; 
	}
	return value; 
}

void Foam::FoamYadeMPI::buildCellPartList(std::shared_ptr<YadeProc>& yProc) {
//  build contribution from particles to the grid from a givebn yade Proc. 
	if (! yProc->foundParticles.size() ) { return; }
	calcInterpWeightGaussian(yProc->foundParticles);
	for (auto& prt : yProc->foundParticles) {
		if (!prt->interpCellWeight.size()) {return; }
		for (unsigned int i =0 ; i != prt->interpCellWeight.size(); ++i){
			const label& cellId = prt->interpCellWeight[i].first; 
			const double& weight = prt->interpCellWeight[i].second; 
			const double& pVol = prt->vol; 
			if ((yProc->pVolContrib.size()==0) || (yProc->uParticleContrib.size()==0)){
				yProc->pVolContrib.push_back(std::make_pair(cellId, pVol*weight)); 
				yProc->uParticleContrib.push_back(std::make_pair(cellId, weight*prt->linearVelocity)); 
				
			} else {
				int c=0; 
				for (unsigned i =0; i != yProc->pVolContrib.size(); ++i) {
					if (yProc->pVolContrib[i].first == cellId && yProc->uParticleContrib[i].first == cellId) {
						yProc->pVolContrib[i].second += (pVol*weight); 
						yProc->uParticleContrib[i].second += (prt->linearVelocity*weight); 
						c +=1; 
					}  
				}
				if (!c) {
						yProc->pVolContrib.push_back(std::make_pair(cellId, pVol*weight)); 
						yProc->uParticleContrib.push_back(std::make_pair(cellId, weight*prt->linearVelocity)); 
					}
			}
		}
	}
}


void Foam::FoamYadeMPI::calcInterpWeightGaussian(std::vector<std::shared_ptr<YadeParticle> >& yParticles ){
  
	// get a list of particles and calculate the gaussian weight 
	if (! yParticles.size()) {return; }
	for (auto & prt : yParticles ){
		prt->interpCellWeight.clear();
		if (! prt->cellIds.size()) {return;  }  
		for (unsigned  i = 0; i !=  prt->cellIds.size(); ++ i){
			double distsq = 0.0; 
			const double& ds1 = mesh.C()[prt->cellIds[i]].x() - prt-> pos.x();
			const double& ds2 = mesh.C()[prt->cellIds[i]].y() - prt-> pos.y();
			const double& ds3 = mesh.C()[prt->cellIds[i]].z() - prt-> pos.z();
			distsq = (ds1*ds1)+ (ds2*ds2) + (ds3*ds3); 
			double weight = exp(-distsq/(2*std::pow(sigmaInterp, 2)))*interpRangeCu*sigmaPi;
			prt -> interpCellWeight.push_back(std::make_pair(prt->cellIds[i], weight)); 
		}
	}  
}

void Foam::FoamYadeMPI::setCellVolFraction(std::shared_ptr<YadeProc>& yProc){
	// set the vol fraction and particle velocity to the grid from the pVolContrib and uParticleContrib arrays from each yade proc. 
	if (!yProc->pVolContrib.size() || !yProc->uParticleContrib.size()) {return; }
	for (unsigned i=0; i != yProc->pVolContrib.size(); ++i){
		const label& id = yProc->pVolContrib[i].first; 
		const double& pvolC = 1.0-(yProc->pVolContrib[i].second/mesh.V()[id]); 
		alpha[id] = ((pvolC > 0) ? pvolC : 1e-08);
		uParticle[id] = yProc->uParticleContrib[i].second/mesh.V()[id]; 
	}
	
}


void Foam::FoamYadeMPI::calcHydroForce(std::shared_ptr<YadeProc>& yProc){
	// calculate hydrodyamic forces, source terms from each yade proc. 
	if (!yProc->foundParticles.size()) {return; }
	for (auto & prt : yProc->foundParticles){
		initParticleForce(prt); 
		if (isGaussianInterp) {
			hydroDragForce(prt); 
			buoyancyForce(prt); 
		} else {
			stokesDragForce(prt); 
		}
		
	}
}

/* Forces */ 

void Foam::FoamYadeMPI::initParticleForce(std::shared_ptr<YadeParticle>& prt){

	prt->hydroForce.x() = small; 
	prt->hydroForce.y() = small; 
	prt->hydroForce.z() = small; 
	
	prt->hydroTorque.x() = small; 
	prt->hydroTorque.y() = small; 
	prt->hydroTorque.z() = small; 
}


void Foam::FoamYadeMPI::hydroDragForce(std::shared_ptr<YadeParticle>&  prt){
	
	if (!prt->interpCellWeight.size()) {return; } 
	vector uf(0,0,0); 
	double alpha_p = 0.0; 
	double pv = 0.0; 
	
	// get velocities and interpolate to particle location.. 
	
	for (unsigned i =0 ; i != prt->interpCellWeight.size(); ++i){
		const int& cellId = prt->interpCellWeight[i].first; 
		const double& wt  = prt->interpCellWeight[i].second; 
		uf = uf + (U[cellId]*wt); 
		alpha_p = alpha_p + ((1.0-alpha[cellId])*wt); 
		pv = pv + (prt->vol*wt); 
	}
	
	alpha_p = alpha_p/(prt->interpCellWeight.size());  
	const double& alpha_f = 1 - alpha_p; 
	
	// force calculation, 
	
	const vector& uRelVel = uf - prt->linearVelocity; 
	const double&  magUR = mag(uRelVel); 
	const double& Re = small + (magUR*prt->dia)/nu; 
	const double& pw = std::pow(Re, 0.687);
	const scalar& Cd = (24/Re)*(1+(0.15*pw));
	const double& coeff = (0.75*Cd*rhoF*magUR)*(1/prt->dia)*(std::pow(alpha_f, -2.65));
	const vector& f = (pv)*coeff*uRelVel*alpha_f;
	prt->hydroForce += f; 
	
	// distribute the drag term 
	
	for (unsigned i =0; i != prt->interpCellWeight.size(); ++i){
		const double& wt =  prt->interpCellWeight[i].second; 
		const int& id = prt->interpCellWeight[i].first; 
		const double& value = (-coeff*wt*(1-alpha[id])); 
		uSourceDrag[id] = uSourceDrag[id] + (value/rhoF); 
		uSource[id] = uSource[id] + (value*uParticle[id]/rhoF); 
	}
}


void Foam::FoamYadeMPI::buoyancyForce(std::shared_ptr<YadeParticle>& prt) {
  
	vector bforce(0,0,0); double pv = 0.0; const vector& gvt = g[0]; 
	for (unsigned i =0; i != prt->interpCellWeight.size(); ++i){
		const double& wt = prt->interpCellWeight[i].second; 
		pv += (prt->vol*wt); 
	}
	bforce = (rhoP-rhoF)*pv*gvt; 
	prt->hydroForce += bforce; 
	
	for (unsigned i =0; i != prt->interpCellWeight.size(); ++i){
		const double& wt = prt->interpCellWeight[i].first; 
		const int& id = prt->interpCellWeight[i].second; 
		const double& ooCellVol = 1./(mesh.V()[id]*rhoF); 
		uSource[id] = uSource[id] + (-bforce*wt*ooCellVol); 
	}
}

void Foam::FoamYadeMPI::addedMassForce(std::shared_ptr<YadeParticle>& prt){
	
	vector ddtUf(0,0,0); double pv = 0.0; 
	
	for (unsigned i=0; i != prt->interpCellWeight.size(); ++i) {
		const double& wt = prt->interpCellWeight[i].second; 
		const int& id = prt->interpCellWeight[i].first; 
		pv += (prt->vol*wt); 
		ddtUf = ddtUf + (ddtU[id]*wt);
	} 
	
	const vector& f = pv*(ddtUf - (prt->linearVelocity/deltaT))*rhoP; 
	prt -> hydroForce += f; 
	
	for (unsigned int i=0; i != prt->interpCellWeight.size(); ++i) {
		const double& wt = prt->interpCellWeight[i].second; 
		const int& id = prt-> interpCellWeight[i].first; 
		const double& ooCellVol = 1./(mesh.V()[id]*rhoF); 
		uSource[id] = uSource[id]+ (-f*wt*ooCellVol); 
	}
	
}

void Foam::FoamYadeMPI::archimedesForce(std::shared_ptr<YadeParticle>& prt) {
	vector divt(0,0,0); vector pg(0,0,0); double pv = 0.0; 
	
	for (unsigned i =0 ; i != prt->interpCellWeight.size(); ++i) {
		const int& cellId = prt->interpCellWeight[i].first; 
		const double& wt = prt->interpCellWeight[i].second; 
		pv +=  (prt->vol*wt); 
		divt = divt + (2.0*nu*divT[cellId]*wt*rhoF);
		pg   =  pg + (-gradP[cellId]*wt);
	}
	
	const vector& f = pv*(divt+pg);
	prt->hydroForce += f; 
	
	for (unsigned i = 0; i != prt->interpCellWeight.size(); ++i){
		const int& cellId = prt->interpCellWeight[i].first; 
		const double& wt = prt->interpCellWeight[i].second; 
		const double& ooCellVol = 1./(mesh.V()[cellId]); 
		uSource[cellId] = uSource[cellId] + (-f*wt*ooCellVol); 
	}
}

void Foam::FoamYadeMPI::stokesDragForce(std::shared_ptr<YadeParticle>& prt) {
  
	autoPtr<interpolation<vector> > interpVel = interpolation<vector>::New("cell", U); 
	const vector& uFluid = interpVel->interpolate(prt->pos, prt->inCell);
	const double& coeff = 3*M_PI*(prt->dia)*nu*rhoF; 
	const double& ooCellVol = 1./(mesh.V()[prt->inCell]*rhoF); 
	prt->hydroForce = coeff*(uFluid-prt->linearVelocity); 
	uSource[prt->inCell]  += (-1*ooCellVol*prt->hydroForce); 
	
}

void Foam::FoamYadeMPI::stokesDragTorque(std::shared_ptr<YadeParticle>& prt) {
  
	autoPtr<interpolation<tensor> > interpGradV = interpolation<tensor>::New("cell", vGrad);  
	const tensor& vGradPt = interpGradV->interpolate(prt->pos, prt->inCell); 
	scalar s1 = vGradPt.zy() - vGradPt.yz(); scalar s2 = vGradPt.zx() - vGradPt.xz(); scalar s3 = vGradPt.yx() - vGradPt.xy(); 
	vector wfluid(s1,s2,s3); 
	prt->hydroTorque = M_PI*(pow(prt->dia, 3))*(wfluid-prt->rotationalVelocity)*nu; 
}

void Foam::FoamYadeMPI::updateSources(std::shared_ptr<YadeProc>& yProc) {
  
	for (unsigned i =0; i != yProc->pVolContrib.size(); ++i) {
		const label& cellId = yProc->pVolContrib[i].first; 
		const double& pvolC = yProc->pVolContrib[i].second; 
		const double& alpha_p = pvolC/mesh.V()[cellId]; 
		uSourceDrag[cellId] = alpha_p*uSourceDrag[cellId]; 
		uSource[cellId] = alpha_p*uSource[cellId]; 
	}
}


void Foam::FoamYadeMPI::calcHydroTorque(std::shared_ptr<YadeProc>& yProc){
	for (auto & prt : yProc->foundParticles) {
		if (isGaussianInterp) {
			scalar s1=0.0; scalar s2=0.0; scalar s3=0.0; 
			for (unsigned i =0; i != prt->interpCellWeight.size(); ++i) {
				const int& cellId = prt->interpCellWeight[i].first; 
				const double& wt = prt->interpCellWeight[i].second; 
				s1 += ((vGrad[cellId].yz() - vGrad[cellId].zy())*wt); 
				s2 += ((vGrad[cellId].zx() - vGrad[cellId].xz())*wt);
				s3 += ((vGrad[cellId].yx() - vGrad[cellId].xy())*wt);
				//s3 += (vGrad[cellId].yx() - vGrad[cellId].xy());
			}
			const vector wfluid(s1,s2,s3); 
			prt->hydroTorque += M_PI*(pow(prt->dia, 3))*(wfluid-prt->rotationalVelocity)*nu*rhoF; 
		} else {
			stokesDragTorque(prt); 
		}
	}
}

// void Foam::FoamYadeMPI::

void Foam::FoamYadeMPI::sendHydroForceYadeMPI(){
	// set forces in hydroForceBuff 
	for (auto& yProc : inCommProcs){
		for (auto & prt : yProc->foundParticles){
			//forces
			yProc->hydroForceBuff[6*(prt->indx)] = prt->hydroForce.x(); 
			yProc->hydroForceBuff[6*(prt->indx+1)] = prt->hydroForce.y();
			yProc->hydroForceBuff[6*(prt->indx+2)] = prt->hydroForce.z();
			//torques
			yProc->hydroForceBuff[6*(prt->indx+3)] = prt->hydroTorque.x();
			yProc->hydroForceBuff[6*(prt->indx+4)] = prt->hydroTorque.y();
			yProc->hydroForceBuff[6*(prt->indx+5)] = prt->hydroTorque.z();
		}  
	}
	// send 
	for (auto & yProc : inCommProcs){
		std::vector<double>& fBuff = yProc->hydroForceBuff; 
		MPI_Send(&fBuff.front(), int(fBuff.size()), MPI_DOUBLE, yProc->yRank, TAG_FORCE, MPI_COMM_WORLD); 
	}
	
}



void Foam::FoamYadeMPI::clearInCommProcs(){
	for (const auto& yProc : inCommProcs){
		yProc->foundParticles.clear(); 
		yProc->numParticles = 0; 
		yProc->pVolContrib.clear(); 
		yProc->uParticleContrib.clear(); 
		yProc->foundParticles.clear(); 
		yProc->hydroForceBuff.clear(); 
		yProc->foundBuff.clear(); 
		yProc->particleDataBuff.clear(); 
		yProc->inComm = false; 
		std::fill(yProc->numParticlesProc.begin(), yProc->numParticlesProc.end(), -1); 
	}
	inCommProcs.clear(); 
	
	for (auto& yProc : yadeProcs){
		std::fill(yProc.numParticlesProc.begin(), yProc.numParticlesProc.end(), -1); 
		yProc.numParticles = 0; 
	}
}

void Foam::FoamYadeMPI::terminateRun(){
	return;
}

void Foam::FoamYadeMPI::calcHydroTimeScale(){
	return; 
}

void Foam::FoamYadeMPI::sendHydroTimeScale(std::shared_ptr<YadeProc>& ){
	// do an all reduce from each yProc. 
	return; 
}

void Foam::FoamYadeMPI::exchangeDT(){
	if (localRank == 0) {
		MPI_Send(&deltaT, 1, MPI_DOUBLE, 0, TAG_FLUID_DT, MPI_COMM_WORLD); 
	}
	
	if (localRank == 0) {
		MPI_Status status; 
		MPI_Recv(&yadeDT, 1, MPI_DOUBLE, 0, TAG_YADE_DT, MPI_COMM_WORLD, &status);  
	}
	// broadcast recvd yadeDt from localRank = 0. 
	MPI_Bcast(&yadeDT,1, MPI_DOUBLE, 0, PstreamGlobals::MPI_COMM_FOAM); 
}


void Foam::FoamYadeMPI::setSourceZero(){
	if(!inCommProcs.size()){return; }
	for (auto yProc : inCommProcs){
		for (auto& pt : yProc->pVolContrib){
			const auto& cellid = pt.first; 
			uSource[cellid].x() = small;
			uSource[cellid].y() = small; 
			uSource[cellid].z() = small; 
			if (isGaussianInterp){
				alpha[cellid] = 1.0; 
				uSourceDrag[cellid] = small; 
				uParticle[cellid].x() = small; 
				uParticle[cellid].y() = small; 
				uParticle[cellid].z() = small; 
			}
			
		}
		yProc->pVolContrib.clear(); 
		yProc->uParticleContrib.clear(); 
	}
	clearInCommProcs();
}


/*main driver*/ 
void Foam::FoamYadeMPI::setParticleAction(double dt) {
	
	deltaT = dt; 
	if (!rankSizeSet) {getRankSize();}
	
	recvYadeIntrs(); 
	locateAllParticles();
	
	if (isGaussianInterp){
		if (!inCommProcs.size()){return; }
		for (auto& yProc : inCommProcs){
			buildCellPartList(yProc);
			setCellVolFraction(yProc);
			calcHydroForce(yProc);
			calcHydroTorque(yProc);
			updateSources(yProc); 
			
		}
	} else {
		for (auto yProc : inCommProcs){
			calcHydroForce(yProc);
			calcHydroTorque(yProc);
		}
	}
	
	sendHydroForceYadeMPI(); 
	exchangeDT(); 
}
