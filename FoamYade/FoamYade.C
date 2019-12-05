/*  (c) 2019  Deepak Kunhappan  deepak.kn1990@gmail.com, deepak.kunhappan@3sr-grenoble.fr */ 

#include "FoamYade.H"
#include "PstreamGlobals.H" 
#include <mpi.h> 
#include <cmath> 
#include <string>

void Foam::FoamYade::setScalarProperties(scalar _rhoP, scalar _rhoF, scalar _nu){
	rhoP = _rhoP; rhoF = _rhoF; nu = _nu; 
}

void Foam::FoamYade::printMsg(const std::string& msg){
	std::cout << " Rank : " << worldRank << "  " << msg << std::endl;   
}


void Foam::FoamYade::getRankSize(){
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
		
		if (commSzDff == 1) {serialYade = true; }
		
		mshTree.build_tree();  // build tree. 
		
		if (!serialYade){
			// alloc vector of yadeProcs, yade Master does not have participate in communications. 
			yadeProcs.resize(commSzDff-1);
			// we do not receive from yade Master, 
			for (unsigned int i=0; i != yadeProcs.size(); ++i) {
				yadeProcs[i].yRank = (int)i+1; 
				yadeProcs[i].numParticlesProc.resize(localCommSize); 
				std::fill(yadeProcs[i].numParticlesProc.begin(), yadeProcs[i].numParticlesProc.end(), -1); 
			}
			sendMeshBbox();
		}
		else {
			YadeProc yProc; 
			yProc.yRank = 0; 
			inCommProcs.push_back(std::make_shared<YadeProc>(yProc)); 
		}
		initFields(); 
	}
}


void Foam::FoamYade::initFields(){

	forAll(uSource,cellI){
		uSource[cellI] = vector(0.0,0.0,0.0); 
		if (isGaussianInterp){
			uParticle[cellI] = vector(0.0,0.0,0.0); 
			alpha[cellI] = 1.0; 
			uSourceDrag[cellI] = 0.0; 
		}
		
	}
	
	alpha = 1.0; 
	interpRange = std::pow(mesh.V()[0], 1.0/3.0);
	if (isGaussianInterp) {interpRange = 3*interpRange; }
	sigmaInterp = interpRange*0.42460; // interp_range/(2sqrt(2ln(2))) filter width half maximum;
	interpRangeCu = std::pow(interpRange, 3.0); 
	sigmaPi = 1.0/(std::pow(2*M_PI*sigmaInterp*sigmaInterp, 1.5));   
}



void Foam::FoamYade::sendMeshBbox(){
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
}


void Foam::FoamYade::recvYadeIntrs(){
  
  
// 	for (int rnk = 0; rnk != commSzDff-1; ++rnk){
// 		MPI_Status status; 
// 		MPI_Recv(&yadeProcs[rnk].numParticlesProc.front(), localCommSize, MPI_INT, rnk+1, TAG_SZ_BUFF, MPI_COMM_WORLD, &status); 
// 	}
  
	for (auto& yProc : yadeProcs){
		MPI_Status status; 
		MPI_Recv(&yProc.numParticlesProc.front(), localCommSize, MPI_INT, yProc.yRank, TAG_SZ_BUFF, MPI_COMM_WORLD, &status); 
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
			
			std::fill(yProc.foundBuff.begin(),yProc.foundBuff.end(), -1); 
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


// serial case
void Foam::FoamYade::allocArrays(int sz, const std::shared_ptr<YadeProc>& yProc){
	yProc->particleDataBuff.clear(); yProc->hydroForceBuff.clear(); 
	if (!fibreCpl) {
		yProc->particleDataBuff.resize(sz*10); 
	} else {
		yProc->particleDataBuff.resize(sz*15); 
	}
	yProc->hydroForceBuff.resize(sz*6);
	yProc->foundBuff.resize(sz);
}



void Foam::FoamYade::locateAllParticles(){
	if (serialYade) {
	  
		  const auto& yProc = inCommProcs[0]; 
		  MPI_Bcast(&yProc->numParticles, 1, MPI_INT, 0, MPI_COMM_WORLD); 
		  // alloc arrys 
		  allocArrays(yProc->numParticles, yProc); 
		  //get particleDataBuff   
		  std::vector<double>& dBuff = yProc->particleDataBuff; 
		  MPI_Bcast(&dBuff.front(), (int) dBuff.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD); 
		  std::fill(yProc->foundBuff.begin(),yProc->foundBuff.end(), -1);
		  sendRanks.resize(yProc->numParticles); std::fill(sendRanks.begin(), sendRanks.end(), -1); 
	}
	
	for (const auto& yProc : inCommProcs){
		for (int np = 0; np != yProc->numParticles; ++np) {
			vector pos(0,0,0); 
			if (!fibreCpl){
				pos.x() = yProc->particleDataBuff[np*10]; 
				pos.y() = yProc->particleDataBuff[np*10+1]; 
				pos.z() = yProc->particleDataBuff[np*10+2]; 
				
			} else {
				pos.x() = yProc->particleDataBuff[np*15];
				pos.y() = yProc->particleDataBuff[np*15+1];
				pos.z() = yProc->particleDataBuff[np*15+2]; 
			}
			
			std::vector<int> cellIds = locatePt(pos); 
			
			int found = 0; 
			
			if (cellIds.size() && cellIds[0] > -1 ){
				std::shared_ptr<YadeParticle> yParticle = std::make_shared<YadeParticle>(); 
				yProc->inComm = true; 
				yParticle->pos = pos; 
				yParticle->cellIds = cellIds; yParticle->inCell = cellIds[0]; 
				yParticle->indx = np; // keep track of 'found particle' index 
				
				yParticle->linearVelocity.x() = yProc->particleDataBuff[np*10+3]; 
				yParticle->linearVelocity.y() = yProc->particleDataBuff[np*10+4]; 
				yParticle->linearVelocity.z() = yProc->particleDataBuff[np*10+5]; 
				
				yParticle->rotationalVelocity.x() = yProc->particleDataBuff[np*10+6]; 
				yParticle->rotationalVelocity.y() = yProc->particleDataBuff[np*10+7]; 
				yParticle->rotationalVelocity.z() = yProc->particleDataBuff[np*10+8]; 
				
				yParticle->dia = 2*yProc->particleDataBuff[np*10+9]; 
				yParticle->calcPartVol(); 
				
				yProc->foundBuff[np] = 1; 
				found = worldRank; 
				yProc->foundParticles.push_back(yParticle); 
			}
			
			if (serialYade) {
				MPI_Allreduce(&found,&sendRanks[np],1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
				if (localRank == 0 && sendRanks[np] < 0) {
					Info << "particle indx = " << np << " pos = " << pos <<  "  was not found in FOAM" << endl; 
				}
			}
		}
		yProc->particleDataBuff.clear(); 
	}
	
	if (!serialYade){
	// send the 'foundBuff' of intersecting procs to Yade. 
		for (const auto& yProc : inCommProcs) {
			std::vector<int>& fBuff = yProc->foundBuff; 
			int sz =  int (fBuff.size()); 
			MPI_Send(&fBuff.front(), sz, MPI_INT, yProc->yRank, TAG_SEARCH_RES, MPI_COMM_WORLD);  
		}
	}
	
}

std::vector<int> Foam::FoamYade::locatePt(const vector& pt){
	if (!isGaussianInterp) {
		std::vector<int> cellId; 
		int inCell = mesh.findCell(pt); 
		if (inCell > -1) cellId.push_back(inCell); 
		return cellId; 
	}
	else {
		return mshTree.nnearestCellsRange(pt, interpRange, true); }
	
}


void Foam::FoamYade::buildCellPartList(YadeProc* yProc) {
//  build contribution from particles to the grid from a givebn yade Proc. 
	if (! yProc->foundParticles.size() ) { return; }
	calcInterpWeightGaussian(yProc->foundParticles);
	for (auto& prt : yProc->foundParticles) {
		for (unsigned int i =0 ; i != prt->interpCellWeight.size(); ++i){
			const label& cellId = prt->interpCellWeight[i].first; 
			const double& weight = prt->interpCellWeight[i].second; 
			const double& pVol = prt->vol; 
			if ((yProc->pVolContrib.size()==0) || (yProc->uParticleContrib.size()==0)){
				yProc->pVolContrib.push_back(std::make_pair(cellId, pVol*weight)); 
				yProc->uParticleContrib.push_back(std::make_pair(cellId, weight*prt->linearVelocity*pVol)); 
				
			} else {
				int c=0; 
				for (unsigned i =0; i != yProc->pVolContrib.size(); ++i) {
					if (yProc->pVolContrib[i].first == cellId && yProc->uParticleContrib[i].first == cellId) {
						yProc->pVolContrib[i].second += (pVol*weight); 
						yProc->uParticleContrib[i].second += (prt->linearVelocity*weight*pVol); 
						c +=1; 
					}  
				}
				if (c==0) {
						yProc->pVolContrib.push_back(std::make_pair(cellId, pVol*weight)); 
						yProc->uParticleContrib.push_back(std::make_pair(cellId, weight*prt->linearVelocity*pVol)); 
					}
			}
		}
	}
}


void Foam::FoamYade::calcInterpWeightGaussian(std::vector<std::shared_ptr<YadeParticle> >& yParticles ){
  
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

void Foam::FoamYade::setCellVolFraction(YadeProc* yProc){
	// set the vol fraction and particle velocity to the grid from the pVolContrib and uParticleContrib arrays from each yade proc. 
	if (!yProc->pVolContrib.size() || !yProc->uParticleContrib.size()) {return; }
	for (unsigned i=0; i != yProc->pVolContrib.size(); ++i){
		const label& id = yProc->pVolContrib[i].first; 
		const double& pvolC = 1.0-(yProc->pVolContrib[i].second/mesh.V()[id]); 
		alpha[id] = ((pvolC > 0) ? pvolC : 1e-08);
		uParticle[id] = yProc->uParticleContrib[i].second/mesh.V()[id]; 
	}
	
}


void Foam::FoamYade::calcHydroForce(YadeProc* yProc){
	// calculate hydrodyamic forces, source terms from each yade proc. 
	if (!yProc->foundParticles.size()) {return; }
	for (const auto & prt : yProc->foundParticles){
		initParticleForce(prt.get()); 
		if (isGaussianInterp) {
			hydroDragForce(prt.get()); 
// 			buoyancyForce(prt); 
		} else {
			stokesDragForce(prt.get()); 
		}
		
	}
}

/* Forces */ 

void Foam::FoamYade::initParticleForce(YadeParticle* prt){
	prt->hydroForce = vector(0,0,0); 
	prt->hydroTorque = vector(0,0,0); 
}


void Foam::FoamYade::hydroDragForce(YadeParticle*  prt){
	
	if (!prt->interpCellWeight.size()) {return; } 
	vector uf(0,0,0); 
	double alpha_p = 0.0; 
	double pv = 0.0; 
	
	// get velocities and interpolate to particle location.. 
	
	for (unsigned i =0 ; i != prt->interpCellWeight.size(); ++i){
		const int& cellId = prt->interpCellWeight[i].first; 
		const double& wt  = prt->interpCellWeight[i].second; 
		uf = uf + (U[cellId]*wt); 
		alpha_p = alpha_p + ((1-alpha[cellId])*wt); 
		pv = pv + (prt->vol*wt); 
	}
	const double& alpha_f = 1 - alpha_p; 
	// force calculation, 
	const vector& uRelVel = uf - prt->linearVelocity; 
	const double&  magUR = mag(uRelVel); 
	const double& Re = small + ((magUR*prt->dia)/nu); 
	vector f; double coeff; 
	
		const double& pw = std::pow(Re, 0.687);
		const scalar& Cd = (24/Re)*(1+(0.15*pw));
		coeff = (0.75*Cd*rhoF*magUR)*(1/prt->dia)*(std::pow(alpha_f, -2.65));
		
		//(0.75*Cd*rhoF*magUR)*(1/prt->dia)*(std::pow(alpha_f, -2.65));
	
	f =  pv*coeff*uRelVel*alpha_f;
	prt->hydroForce += f;
	
	// distribute the drag term 
	
	for (unsigned i =0; i != prt->interpCellWeight.size(); ++i){
		const double& wt =  prt->interpCellWeight[i].second; 
		const int& id = prt->interpCellWeight[i].first; 
		const double& value = (-coeff*wt*(1-alpha[id])); 
		uSourceDrag[id] = uSourceDrag[id] + (value/rhoF); 
		uSource[id] = uSource[id] + ((value*uParticle[id])/rhoF); 
	}
}


void Foam::FoamYade::buoyancyForce(YadeParticle* prt) {
  
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

void Foam::FoamYade::addedMassForce(YadeParticle* prt){
	
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

void Foam::FoamYade::archimedesForce(YadeParticle* prt) {
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

void Foam::FoamYade::stokesDragForce(YadeParticle* prt) {
	autoPtr<interpolation<vector> > interpVel = interpolation<vector>::New("cell", U); 
	const vector& uFluid = interpVel->interpolate(prt->pos, prt->inCell);
	const double& coeff = 3*M_PI*(prt->dia)*nu*rhoF; 
	const double& ooCellVol = 1./(mesh.V()[prt->inCell]*rhoF); 
	prt->hydroForce = coeff*(uFluid-prt->linearVelocity); 
	uSource[prt->inCell]  += (-1*ooCellVol*prt->hydroForce); 
}

void Foam::FoamYade::stokesDragTorque(YadeParticle* prt) {
  
	autoPtr<interpolation<tensor> > interpGradV = interpolation<tensor>::New("cell", vGrad);  
	const tensor& vGradPt = interpGradV->interpolate(prt->pos, prt->inCell); 
	scalar s1 = vGradPt.zy() - vGradPt.yz(); scalar s2 = vGradPt.zx() - vGradPt.xz(); scalar s3 = vGradPt.yx() - vGradPt.xy(); 
	vector wfluid(s1,s2,s3); 
	prt->hydroTorque = M_PI*(pow(prt->dia, 3))*(wfluid-prt->rotationalVelocity)*nu; 
}

void Foam::FoamYade::updateSources(YadeProc* yProc) {
  
	for (unsigned i =0; i != yProc->pVolContrib.size(); ++i) {
		const label& cellId = yProc->pVolContrib[i].first; 
		const double& pvolC = yProc->pVolContrib[i].second; 
		const double& alpha_p = pvolC/mesh.V()[cellId]; 
		uSourceDrag[cellId] = alpha_p*uSourceDrag[cellId]; 
		uSource[cellId] = alpha_p*uSource[cellId]; 
	}
}


void Foam::FoamYade::calcHydroTorque(YadeProc* yProc){
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
			stokesDragTorque(prt.get()); 
		}
	}
}

// void Foam::FoamYade::

void Foam::FoamYade::sendHydroForceYadeMPI(){
	// set forces in hydroForceBuff 
	for (auto& yProc : inCommProcs){
		for (auto & prt : yProc->foundParticles){
			//forces
			yProc->hydroForceBuff[6*prt->indx]   = prt->hydroForce.x(); 
			yProc->hydroForceBuff[6*prt->indx+1] = prt->hydroForce.y();
			yProc->hydroForceBuff[6*prt->indx+2] = prt->hydroForce.z();
			//torques
			yProc->hydroForceBuff[6*prt->indx+3] = prt->hydroTorque.x();
			yProc->hydroForceBuff[6*prt->indx+4] = prt->hydroTorque.y();
			yProc->hydroForceBuff[6*prt->indx+5] = prt->hydroTorque.z();
		}  
	}
	
	if (!serialYade){
	// send 
		for (auto & yProc : inCommProcs){
			std::vector<double>& fBuff = yProc->hydroForceBuff; 
			MPI_Send(&fBuff.front(), int(fBuff.size()), MPI_DOUBLE, yProc->yRank, TAG_FORCE, MPI_COMM_WORLD); 
		} 
	} else {
	//all reduce
		if (isGaussianInterp)
			for (int np = 0; np != inCommProcs[0]->numParticles; ++np) {
				for (int j = 0; j != 6; ++j) {
				double dummy = 0.0; 
				MPI_Allreduce(&inCommProcs[0]->hydroForceBuff[6*np+j],&dummy ,1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			}
		}
		else {
			const auto& yProc = inCommProcs[0]; 
			for (int np = 0; np != yProc->numParticles; ++np){
				if (sendRanks[np] == worldRank) {
					std::vector<double> fBuff = {yProc->hydroForceBuff[6*np],
								    yProc->hydroForceBuff[6*np+1],
								    yProc->hydroForceBuff[6*np+2],
								    yProc->hydroForceBuff[6*np+3],
								    yProc->hydroForceBuff[6*np+4],
								    yProc->hydroForceBuff[6*np+5] }; 
					
					
					MPI_Send(&fBuff.front(), 6, MPI_DOUBLE, 0, TAG_FORCE, MPI_COMM_WORLD); 
				}
			}
		}
	}
	
}

void Foam::FoamYade::exchangeDT(){
	if (localRank == 0) {
		MPI_Send(&deltaT, 1, MPI_DOUBLE, 0, TAG_FLUID_DT, MPI_COMM_WORLD); 
	}
	if (!serialYade){
		if (localRank == 0) {
			MPI_Status status; 
			MPI_Recv(&yadeDT, 1, MPI_DOUBLE, 0, TAG_YADE_DT, MPI_COMM_WORLD, &status);  
		}
		// broadcast recvd yadeDt from localRank = 0. 
		MPI_Bcast(&yadeDT,1, MPI_DOUBLE, 0, PstreamGlobals::MPI_COMM_FOAM); 
	} else {
		MPI_Bcast(&yadeDT,1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
	}
	
	
}


void Foam::FoamYade::setSourceZero(){
	forAll(uSource, cellI){
		uSource[cellI] = vector(0,0,0); 
		if (isGaussianInterp){
			alpha[cellI] = 1.0; 
			uSourceDrag[cellI] = 0.0;
			uParticle[cellI] = vector(0.0,0.0,0.0); 
		}
	}
	clearInCommProcs();
}

void Foam::FoamYade::clearInCommProcs(){
	if (inCommProcs.size()){
		for (const auto& yp : inCommProcs){
			yp->foundBuff.clear(); 
			yp->foundParticles.clear(); 
			yp->hydroForceBuff.clear(); 
			yp->pVolContrib.clear(); 
			yp->uParticleContrib.clear(); 
		}
		if (!serialYade) inCommProcs.clear(); 
	}
  
}

//TODO
void Foam::FoamYade::calcHydroTimeScale(){
	return; 
}

void Foam::FoamYade::sendHydroTimeScale(YadeProc* yProc ){
	// do an all reduce from each yProc. 
	return; 
}




void Foam::FoamYade::finalizeRun(){
	int value = -1; 
	MPI_Bcast(&value, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (value == 10) MPI_Finalize(); 
}




/*main driver*/ 
void Foam::FoamYade::setParticleAction(double dt) {
	
	deltaT = dt; 
	if (!rankSizeSet) {getRankSize();}
	if (!serialYade) recvYadeIntrs(); 
	locateAllParticles();
	
	if (isGaussianInterp){
		if (inCommProcs.size()){
			for (const auto& yProc : inCommProcs){
				buildCellPartList(yProc.get());
				setCellVolFraction(yProc.get());
				calcHydroForce(yProc.get());
				calcHydroTorque(yProc.get());
			}
		}
	} else {
		if (inCommProcs.size()){
			for (const auto yProc : inCommProcs){
				calcHydroForce(yProc.get());
				calcHydroTorque(yProc.get());
			}
		}
	}
	
	sendHydroForceYadeMPI(); 
	exchangeDT(); 
}
