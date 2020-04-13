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
		// alloc vector of yadeProcs, yade Master does not have participate in communications. 
		
		if (!serialYade) {
			yadeProcs.resize(commSzDff-1);
			// we do not receive from yade Master, 
			for (unsigned int i=0; i != yadeProcs.size(); ++i) {
				yadeProcs[i].yRank = (int)i+1; 
				yadeProcs[i].numParticlesProc.resize(localCommSize); 
				std::fill(yadeProcs[i].numParticlesProc.begin(), yadeProcs[i].numParticlesProc.end(), -1); 
			}
		} else {
			yadeProcs.resize(1); 
			yadeProcs[0].yRank = 0; 
			yadeProcs[0].numParticlesProc.resize(localCommSize); 
			std::fill(yadeProcs[0].numParticlesProc.begin(), yadeProcs[0].numParticlesProc.end(), -1); 
		}
		sendMeshBbox();
		initFields(); 
	}
}


void Foam::FoamYade::initFields(){

	forAll(uSource,cellI){
		uSource[cellI] = vector(0.0,0.0,0.0); 
		if (isGaussianInterp){
			uParticle[cellI] = vector(0.0,0.0,0.0); 
			alpha[cellI] = 1.0; 
			uSourceDrag[cellI] =0.0;
			uCoeff[cellI] = 0.0; 
			uInterp[cellI] = vector(0.0,0.0,0.0); 
		}
		
	}
	alpha = 1.0; 
	interpRange = 1.0*std::pow(mesh.V()[0], 1.0/3.0);
	//sigmaInterp = interpRange*0.42460; // interp_range/(2sqrt(2ln(2))) filter width half maximum;
	//interpRangeCu = std::pow(interpRange, 3.0); 
	sigmaPi = 1.0/(std::pow(M_PI*interpRange*interpRange, 1.5));   
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


void Foam::FoamYade::locateAllParticles(){
	
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
				yParticle->indx = np; 
				
				yProc->foundBuff[np] = 1; 
				
				yProc->foundParticles.push_back(yParticle); 
			} 
		}
		yProc->particleDataBuff.clear(); 
	}
	
	// send the 'foundBuff' of intersecting procs to Yade. 
	for (const auto& yProc : inCommProcs) {
		std::vector<int>& fBuff = yProc->foundBuff; 
		int sz =  int (fBuff.size()); 
		MPI_Send(&fBuff.front(), sz, MPI_INT, yProc->yRank, TAG_SEARCH_RES, MPI_COMM_WORLD);  
	}
	
	//recv info of 'shared particles'. 
	std::vector<int> sharedCount; std::vector<int> rnks; 
	sharedCount.resize(inCommProcs.size()); 

	int n = 0; 
	for (const auto& yProc : inCommProcs) {
		int& ncount  = sharedCount[n]; 
		MPI_Status stat; 
		MPI_Recv(&ncount, 1, MPI_INT, yProc->yRank, TAG_SHARED_ID, MPI_COMM_WORLD, &stat);  
		if (ncount > 0)  rnks.push_back(n);
		++n; 
	}
	

	 std::vector<std::vector<std::vector<int>>> sharedBuff; sharedBuff.resize(inCommProcs.size());
	// prealloc some arrays in sharedBuff; 
	for (unsigned i =0; i != sharedCount.size(); ++i) {
		if (sharedCount[i] > 0) sharedBuff[i].resize(sharedCount[i]); 
	}
	
		
	 
	if (rnks.size()) {
	  
		for (unsigned i = 0; i != rnks.size();  ++i) {
			const int& sender  = inCommProcs[rnks[i]]->yRank; 
			const int& ncount  = sharedCount[rnks[i]]; 
			if (ncount == 0) continue; 
			for (int ii =0; ii != ncount; ++ ii) {
				std::vector<int> buff; 
				MPI_Status status;
				MPI_Probe(sender, TAG_SHARED_ID, MPI_COMM_WORLD, &status);
				int sz; 
				MPI_Get_count(&status, MPI_INT, &sz);
				buff.resize(sz); 
				MPI_Recv(&buff.front(), sz, MPI_INT, sender, TAG_SHARED_ID, MPI_COMM_WORLD, &status);  
				sharedBuff[rnks[i]][ii] = buff; 
			}
		}
// 		// keep indexes of 'shared particle' from each yade proc.
// 		for (unsigned i =0; i != rnks.size(); ++i) {
// 			const std::shared_ptr<YadeProc>& yPrc = inCommProcs[rnks[i]]; 
// 			const int& ncount = sharedCount[rnks[i]]; 
// 			for (int ii=0; ii != ncount; ++ii) {
// 				const auto& buff = sharedBuff[rnks[i]][ii]; 
// 				const int& pIndx = buff[0]; 
// 				const int& shardIndx = getpIndx(yPrc->foundParticles, pIndx); 
// 				if (shardIndx >=  0 ) {
// 					yPrc->foundParticles[shardIndx]->sharedParticle = true ; 
// 					yPrc->foundParticles[shardIndx]->commRanks.resize(buff.size()-1); 
// 					std::copy(buff.begin()+1, buff.end(), yPrc->foundParticles[shardIndx]->commRanks.begin()); 
// 					yPrc->sharedPindxs.push_back(shardIndx);
// 				}
// 				
// 				
// 			}
// 		}
// 		
// // 		//define 'owner' for the shared particles. 
// 		for (const auto& yProc : inCommProcs) {
// 			for (const auto& inx  : yProc->sharedPindxs) {
// 				const auto& prt = yProc->foundParticles[inx];
// 				int numCells = prt->cellIds.size(); 
// 				std::vector<int> otherCellSz; otherCellSz.resize(prt->commRanks.size());
// 				// send and recv the sizes of cellids from each proc. The proc having the largest sz of cellids is the  particle owner. 
// 				// first send.. 
// 				std::vector<MPI_Request> reqs; 
// 				
// 				for (const auto& rnk : prt->commRanks){
// 					MPI_Request req; 
// 					MPI_Isend(&numCells, 1, MPI_INT, rnk, TAG_ID, MPI_COMM_WORLD, &req); 
// 					reqs.push_back(req);
// 				}
// 				// recv . 
// 				int n = 0; 
// 				for (const auto& rnk : prt->commRanks) {
// 					MPI_Status stat; 
// 					MPI_Recv(&otherCellSz[n], 1, MPI_INT, rnk, TAG_ID, MPI_COMM_WORLD, &stat); 
// 					++n; 
// 				}
// 				
// 				// complete MPI_Isend 
// 				for (auto& req : reqs) {
// 					MPI_Status stat; 
// 					MPI_Wait(&req, &stat); 
// 				}
// 				
// 			}
// 		} 
// 		
	}
}


int Foam::FoamYade::getpIndx(const std::vector<std::shared_ptr<YadeParticle>>& particles, const int& indx){
	int inx = -1; 
	int n = 0; 
	for (const auto& prt : particles) {
		if (prt->indx == indx){inx = n; return inx;  }
		++n; 
	}
	return inx; 
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
			setPvolContrib(cellId, prt->vol*weight, prt->linearVelocity*weight*prt->vol); 
		}
	}
	
}

void Foam::FoamYade::setPvolContrib(const int& cellId, const double& pVol, const vector& pVel) {
	// if empty insert contribution, 
	if (! pVolcontrib.size() ) {
			pVolcontrib.insert({cellId, std::make_pair(pVol, pVel)}); 
			cellCount.push_back(std::make_pair(cellId, 1));
			return;  
	}
	std::map<int, std::pair<double, vector> >::iterator it = pVolcontrib.find(cellId);
	if ( it == pVolcontrib.end()) {
		pVolcontrib.insert({cellId, std::make_pair(pVol, pVel)});   
		cellCount.push_back(std::make_pair(cellId,1));
	} else { 
		
		it->second.first += pVol; 
		it->second.second += pVel; 
		incrementCellCount(cellId); 

	}
}

void Foam::FoamYade::incrementCellCount(const int& cellId) {
	for (auto& cellC : cellCount) {
		if (cellC.first == cellId) {cellC.second += 1;} 
	}
}


void Foam::FoamYade::setCellVolFraction(){
	// set the fluid vol fraction and particle velocity to the grid; 
	if (!pVolcontrib.size()) {return; } 
	for (const auto& iv : pVolcontrib) {
		const int& cellId = iv.first; 
		const double& pVol = iv.second.first/mesh.V()[cellId]; 
		const vector& pVel = iv.second.second/mesh.V()[cellId]; 
		const double& vfrc = 1.0 - pVol; 
		alpha[cellId] = (vfrc > 0.20 && vfrc < 1.0 ) ? vfrc : 0.20; 
		uParticle[cellId] = pVel/(1.0-alpha[cellId]); 
		
	}
	
}

void Foam::FoamYade::calcInterpWeightGaussian(std::vector<std::shared_ptr<YadeParticle> >& yParticles ){
  
	// get a list of particles and calculate the gaussian weight 
	if (! yParticles.size()) {return; }
	for (auto & prt : yParticles ){
		
		prt->interpCellWeight.clear();
		if (! prt->cellIds.size()) {return;  }  
		double allwt = 0.0; 
		for (unsigned  i = 0; i !=  prt->cellIds.size(); ++ i){
			double distsq = 0.0; 
			const double& ds1 = mesh.C()[prt->cellIds[i]].x() - prt-> pos.x();
			const double& ds2 = mesh.C()[prt->cellIds[i]].y() - prt-> pos.y();
			const double& ds3 = mesh.C()[prt->cellIds[i]].z() - prt-> pos.z();
			distsq = (ds1*ds1)+ (ds2*ds2) + (ds3*ds3); 
			double weight = exp(-distsq/(2*std::pow(interpRange, 2)))*sigmaPi;
			allwt += weight; 
			prt -> interpCellWeight.push_back(std::make_pair(prt->cellIds[i], weight)); 
		}
		for (auto& it : prt->interpCellWeight){
			it.second = it.second/allwt; 
		}
	}  
}




// void Foam::FoamYade::buildCellPartList(YadeProc* yProc) {
// //  build contribution from particles to the grid from a givebn yade Proc. 
// 	if (! yProc->foundParticles.size() ) { return; }
// 	calcInterpWeightGaussian(yProc->foundParticles);
// 	for (auto& prt : yProc->foundParticles) {
// 		for (unsigned int i =0 ; i != prt->interpCellWeight.size(); ++i){
// 			const label& cellId = prt->interpCellWeight[i].first; 
// 			const double& weight = prt->interpCellWeight[i].second; 
// 			const double& pVol = prt->vol; 
// 			if ((yProc->pVolContrib.size()==0) || (yProc->uParticleContrib.size()==0)){
// 				yProc->pVolContrib.push_back(std::make_pair(cellId, pVol*weight)); 
// 				yProc->uParticleContrib.push_back(std::make_pair(cellId, weight*prt->linearVelocity*pVol)); 
// 				
// 			} else {
// 				int c=0; 
// 				for (unsigned i =0; i != yProc->pVolContrib.size(); ++i) {
// 					if (yProc->pVolContrib[i].first == cellId && yProc->uParticleContrib[i].first == cellId) {
// 						yProc->pVolContrib[i].second += (pVol*weight); 
// 						yProc->uParticleContrib[i].second += (prt->linearVelocity*weight*pVol); 
// 						c +=1; 
// 					}  
// 				}
// 				if (c==0) {
// 						yProc->pVolContrib.push_back(std::make_pair(cellId, pVol*weight)); 
// 						yProc->uParticleContrib.push_back(std::make_pair(cellId, weight*prt->linearVelocity*pVol)); 
// 					}
// 			}
// 		}
// 	}
// }





// void Foam::FoamYade::setCellVolFraction(YadeProc* yProc){
// 	// set the vol fraction and particle velocity to the grid from the pVolContrib and uParticleContrib arrays from each yade proc. 
// 	if (!yProc->pVolContrib.size() || !yProc->uParticleContrib.size()) {return; }
// 	for (unsigned i=0; i != yProc->pVolContrib.size(); ++i){
// 		const label& id = yProc->pVolContrib[i].first; 
// 		const double& pvolC = 1.0-(yProc->pVolContrib[i].second/mesh.V()[id]); 
// 		Pout << "pvolC = " << pvolC << endl; 
// 		alpha[id] = ((pvolC > 0.10) ? pvolC : 0.10);
// 		uParticle[id] = yProc->uParticleContrib[i].second/(mesh.V()[id]*(1-alpha[id])); 
// 	}
// 	
// }


void Foam::FoamYade::calcHydroForce(YadeProc* yProc){
	// calculate hydrodyamic forces, source terms from each yade proc. 
	if (!yProc->foundParticles.size()) {return; }
	for (const auto & prt : yProc->foundParticles){
		initParticleForce(prt.get()); 
		if (isGaussianInterp) {
			hydroDragForce(prt.get()); 
			archimedesForce(prt.get()); 
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
  
	
	if (!prt->interpCellWeight.size()){return; } 
	vector uf(0,0,0); 
	double alpha_f = 0.0; 
	double pv = 0.0; 
	
	for (const auto& idwt : prt->interpCellWeight){
		uf += (U[idwt.first]*idwt.second); 
		alpha_f += (alpha[idwt.first]*idwt.second); 
		pv += (prt->vol*idwt.second); 
	}

	double alpha_p = 1-alpha_f; 
	
	const vector urelvel = (uf-prt->linearVelocity); 
	const double magUR = mag(urelvel); 
	const double Re =  ((magUR*prt->dia)/nu);
	double cd = 0;  
	if (Re > small) cd = Re < 1000 ? (24/(Re))*(1+(0.15*std::pow(Re, 0.687))) : 0.44; 
	double coeff; 
	
	if (alpha_f > 0.8){
		coeff = 0.75*cd*alpha_f*alpha_p*(1/prt->dia)*rhoF*magUR*std::pow(alpha_f, -2.65); 
	} else {
		const double& mu_f = nu*rhoF; 
		const double& alpha_pSq = std::pow(alpha_p, 2); 
		const double& dpSq = std::pow(prt->dia, 2); 
		double cf1 = 150*mu_f*(alpha_pSq)*(1/(alpha_f*dpSq)); 
		double cf2 = 1.75*alpha_p*rhoF*(1/(prt->dia*alpha_f))*magUR; 
		coeff = cf1+cf2; 
	}
	

	vector hf = pv*coeff*urelvel*(1.0/alpha_p); 
	prt->hydroForce += hf; 
	
	for (const auto& idwt : prt->interpCellWeight){
		double ooCellVol = (prt->vol)/(rhoF*mesh.V()[idwt.first]); 
		uInterp[idwt.first] +=  (uf*idwt.second); 
		uSourceDrag[idwt.first] += (-coeff*idwt.second*ooCellVol); 
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
	pv = pv/prt->interpCellWeight.size(); 
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
		pg   =  pg + (gradP[cellId]*wt);
	}
	
	const vector& f = pv*(-pg-divt);
	prt->hydroForce += f; 
	
	for (unsigned i = 0; i != prt->interpCellWeight.size(); ++i){
		const int& cellId = prt->interpCellWeight[i].first; 
		const double& wt = prt->interpCellWeight[i].second; 
		const double& ooCellVol = 1./(mesh.V()[cellId]*rhoF); 
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
	prt->hydroTorque = M_PI*(pow(prt->dia, 3))*(wfluid-prt->rotationalVelocity)*nu*rhoF; 
}

void Foam::FoamYade::updateSources() {

	if (!pVolcontrib.size()) {return; } 
	
	for (const auto& iv : pVolcontrib) {
		const int& cellId = iv.first;
		const double alphaf = alpha[cellId]; const double alphap = 1-alpha[cellId]; 
		Pout << "Ucell = " << U[cellId] << endl; 
		Pout << "uInterp = " << uInterp[cellId] << endl; 
		Pout << "uParticle = " << uParticle[cellId] << endl;  
		uCoeff[cellId] = mag(U[cellId]) > small ? (U[cellId] & uInterp[cellId])/(mag(U[cellId])) : 0.0;  
		uSourceDrag[cellId] = uSourceDrag[cellId]/(alphap*alphaf);
		uCoeff[cellId] = uCoeff[cellId]*uSourceDrag[cellId];  
		Pout << "ucoeff = "  << uCoeff[cellId] << endl; 
		Pout << "uDrag = "  << uSourceDrag[cellId] << endl; 
	}
/*	
	for (const auto& cellC : cellCount) {
		uSourceDrag[cellC.first] = uSourceDrag[cellC.first]/(cellC.second); 
	}*/

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
			//torquesupdate
			yProc->hydroForceBuff[6*prt->indx+3] = prt->hydroTorque.x();
			yProc->hydroForceBuff[6*prt->indx+4] = prt->hydroTorque.y();
			yProc->hydroForceBuff[6*prt->indx+5] = prt->hydroTorque.z();
		}  
	}
	
	for (auto & yProc : inCommProcs){
		std::vector<double>& fBuff = yProc->hydroForceBuff; 
		MPI_Send(&fBuff.front(), int(fBuff.size()), MPI_DOUBLE, yProc->yRank, TAG_FORCE, MPI_COMM_WORLD); 
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
			uInterp[cellI] = vector(0.0,0.0,0.0); 
			uCoeff[cellI] = 0.0; 
		}
		
	}
	pVolcontrib.clear(); 
	cellCount.clear(); 
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
		inCommProcs.clear();  
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
	recvYadeIntrs(); 
	locateAllParticles();
	
	if (isGaussianInterp){
		if (inCommProcs.size()){
			for (const auto& yProc : inCommProcs){
				buildCellPartList(yProc.get());
				//setCellVolFraction(yProc.get());
				//calcHydroForce(yProc.get());
				//updateSources(yProc.get());
 				//calcHydroTorque(yProc.get());
			}
			setCellVolFraction();
			for (const auto& yProc : inCommProcs) {
				calcHydroForce(yProc.get());
			}
			updateSources();
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
