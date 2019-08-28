//Deepak Kn , deepak.kunhappan@3sr-grenoble.fr, deepak.kn1990@gmail.com Feb 2019
// Part of Yade-OpenFOAM coupling.
// TODO:  calculate relaxation time for the particles based on drag coeff.

#include "foamYade.H"
#include <iostream>
#include <memory>
#include <cmath>

/************************************************************************************/

void Foam::foamYade::initAllocArrays(){

  comm.cast_integer_data(yadeProc, numParticles);
  hydroForce.reserve(numParticles*6);
  particleData.reserve(numParticles*10); // particleLinearVel, rotVel, ori, dia;
  particleInProc.reserve(numParticles);
  particleInProc.assign(numParticles,-1);
  hydroForce.assign(numParticles*6,1e-19);
  haveParticle= false;
  recvdParticleData = false;
  mshTree.build_tree();
  interp_range = 2*std::pow(mesh.V()[0], 1.0/3.0); // assuming uniform mesh, change later.
  sigma_interp = interp_range*0.42460; // interp_range/(2sqrt(2ln(2))) filter width half maximum;
  sigma_pi = 1.0/(std::pow(2*M_PI*sigma_interp*sigma_interp, 1.5));
  interp_range_cu = interp_range*interp_range*interp_range; //
  //Initialize alpha
  alpha = 1.0;
  forAll(uSource, cellI){
    uSource[cellI].x() = 1e-15;
    uSource[cellI].y() = 1e-15;
    uSource[cellI].z() = 1e-15;
    uSourceDrag[cellI] = 1e-15;
    uParticle[cellI].x() = 1e-15;
    uParticle[cellI].y() = 1e-15;
    uParticle[cellI].z() = 1e-15; //  no std::fill
  }
}

/************************************************************************************/


 void Foam::foamYade::allocArrays(int sz){
    hydroForce.resize(sz*6);
    particleData.resize(sz*6); 
    particleInProc.resize(sz*6); 
    forAll(uSource, cellI){
    uSource[cellI].x() = 1e-15;
    uSource[cellI].y() = 1e-15;
    uSource[cellI].z() = 1e-15;
    uSourceDrag[cellI] = 1e-15;
    uParticle[cellI].x() = 1e-15;
    uParticle[cellI].y() = 1e-15;
    uParticle[cellI].z() = 1e-15; //  no std::fill
  }
}


/************************************************************************************/

void Foam::foamYade::calcInterpWeightGaussian(std::vector<yadeParticle*>& localParticle) {

  // calculate weights for averaging based on Gaussian distribution.

  if (localParticle.size()==0) {return; }

  for (std::vector<yadeParticle*>::iterator pIter=localParticle.begin(); pIter != localParticle.end(); ++pIter){

    yadeParticle* particle = *pIter;

     if (! particle -> cellIds.size()) {return ; }

     for (unsigned int i=0; i != particle -> cellIds.size(); ++i ) {
       double distsq = 0;
       const double& ds1 = mesh.C()[particle -> cellIds[i]].x() - particle -> pos.x();
       const double& ds2 = mesh.C()[particle -> cellIds[i]].y() - particle -> pos.y();
       const double& ds3 = mesh.C()[particle -> cellIds[i]].z() - particle -> pos.z();
       distsq = (ds1*ds1)+(ds2*ds2)+(ds3*ds3);
       double weight = exp(-distsq/(2*std::pow(sigma_interp, 2)))*interp_range_cu*sigma_pi;
       particle -> interpCellWeight.push_back(std::make_pair(particle-> cellIds[i],weight));
     }
       //sum the weights
//      double wSum = 0.0;
//      for (unsigned int i=0; i != particle -> interpCellWeight.size(); ++i ) {
//        wSum += particle->interpCellWeight[i].second;  }
//       // divide by sum:
//      for (unsigned int i=0; i != particle -> interpCellWeight.size(); ++i ) {
//        particle -> interpCellWeight[i].second = particle->interpCellWeight[i].second/wSum;
//      }
   }
}

/************************************************************************************/

void Foam::foamYade::setScalarProperties(scalar nu_val, scalar pDensity, scalar fDensity) {
  nu = nu_val;
  partDensity = pDensity;
  fluidDensity = fDensity;
}

/************************************************************************************/

void Foam::foamYade::buildCelltoPartList(std::vector<std::pair<label, double> >& pVolContrib, std::vector<std::pair<label, vector> >& uParticleContrib) {

  // we need the weighted particle data on the grid.
  //pVolContrib : vector (c++) consisting of the cell ids and weighted volume (summed)
  //uParticleContrib : vector (c++) consisting of the cell ids and weighted particle velocity (summed)
  //maybe just one list of ids are enough?


  //Return if this proc has no particles.

  if (localParticle.size() ==0 ) {return; }

  else {
    calcInterpWeightGaussian(localParticle); // get weights of the corresponding particle - cell ids

    // loop throught the local particle list .

    for (std::vector<yadeParticle*>::iterator pIter = localParticle.begin(); pIter != localParticle.end(); ++pIter) {

      // loop through the particle's list of cell ids it's located.

      yadeParticle* particle = *pIter;
      for (unsigned int i=0; i != particle->interpCellWeight.size(); ++i ) {

        const label& cellid = particle-> interpCellWeight[i].first;
        const double& weight = particle-> interpCellWeight[i].second;
        const double& pvol =  particle-> vol;

       // if the contribution vectors are empty, fill them.

       if (pVolContrib.size()==0 || uParticleContrib.size() ==0) {
         pVolContrib.push_back(std::make_pair(cellid, pvol*weight));
         uParticleContrib.push_back(std::make_pair(cellid,particle -> linearVelocity*weight*pvol));
       }
        //check if the particle->cellid exists in the pVolContrib.first
        else {
         int c =0;
         for (unsigned int i=0; i != pVolContrib.size(); ++i) {
            if (pVolContrib[i].first == cellid || uParticleContrib[i].first ==cellid) {
              // if found include this contribution.
              pVolContrib[i].second += (pvol*weight);
              uParticleContrib[i].second += (particle->linearVelocity*weight*pvol);
              c +=1;
            }
            // if not previously added, add now.
         } if (c==0){ pVolContrib.push_back(std::make_pair(cellid, pvol*weight));
                      uParticleContrib.push_back(std::make_pair(cellid, particle -> linearVelocity*weight*pvol));
                } }
      }
    }
  }

}

/*********************************************************************************************/

void Foam::foamYade::resetLocalParticleList(std::vector<yadeParticle*>& localParticle)
{
  const std::size_t& vec_sz = localParticle.size();
  if (vec_sz > 0){
    while (!localParticle.empty()) delete localParticle.back(), localParticle.pop_back(); }
  localParticle.clear();
}

/*********************************************************************************************/

void Foam::foamYade::locateAllParticle()
{
       MPI_Bcast(&numParticles, 1, MPI_INT, 0, MPI_COMM_WORLD);    
       allocArrays(numParticles);  
       
       comm.cast_double_array_data(yadeProc, particleData);
       
       
       
       for (int i=0; i != numParticles; ++i)
       {

         yadeParticle*  particle = new yadeParticle();
         //std::unique_ptr<yadeparticle-> particle (new yadeParticle);
         particle->pos.x() = particleData[i*10];
         particle->pos.y()= particleData[i*10+1];
         particle->pos.z() = particleData[i*10+2];
         particle->inProc = -1;

        if (locateParticle(particle)){

           particle-> inProc = comm.rank;
           particle->indx = i;
           localParticle.push_back(particle);
           particle->linearVelocity.x()=particleData[i*10+3];
           particle->linearVelocity.y()= particleData[i*10+4];
           particle->linearVelocity.z() = particleData[i*10+5];
           particle->rotationalVelocity.x()=particleData[i*10+6];
           particle->rotationalVelocity.y()=particleData[i*10+7];
           particle->rotationalVelocity.z()=particleData[i*10+8];
           particle->dia = 2*particleData[i*10+9];
           particle -> calcPartVol(2*particleData[i*10+9]);

        }
        comm.procReduceMaxInt(particle-> inProc, particleInProc[i]);
        if (comm.rank==1 && particleInProc[i] < 0)
           std::cout << "particle->id = " << i  << " " << " pos " << particle -> pos.x() << " " << particle -> pos.y() << " "<< particle -> pos.z() << " "<< "not found " << std::endl;
        if (particle -> inProc != comm.rank)
          delete particle;
  }
}
/************************************************************************************/

bool Foam::foamYade::locateParticle(yadeParticle* aParticle){

    bool value = false;
   if (isGaussianInterp ){
    aParticle ->cellIds = mshTree.nnearestCellsRange(aParticle->pos, interp_range, isGaussianInterp);
    if (aParticle->cellIds.size() > 0 ){value = true; }
   } else {
     label cellid = mesh.findCell(aParticle -> pos);
     if (cellid > -1 ) { aParticle->inCell= cellid; value = true;  }
   }
    return value;
}
/************************************************************************************/

void Foam::foamYade::sendHydroForcePoint() {

  for (int i=0; i != numParticles; ++i)
  {
    if (comm.rank==particleInProc[i])
    {
      std::size_t vec_sz = localParticle.size();
      for (unsigned int ii=0; ii!=vec_sz; ++ii)
      {
         if(localParticle[ii]->indx== i) {
           comm.sendOneDouble(yadeProc,localParticle[ii]->hydroForce.x());
           comm.sendOneDouble(yadeProc,localParticle[ii]->hydroForce.y());
           comm.sendOneDouble(yadeProc,localParticle[ii]->hydroForce.z());
           comm.sendOneDouble(yadeProc,localParticle[ii]->hydroTorque.x());
           comm.sendOneDouble(yadeProc,localParticle[ii]->hydroTorque.y());
           comm.sendOneDouble(yadeProc,localParticle[ii]->hydroTorque.z());
         }
      }
    }
  }
  
  resetLocalParticleList(localParticle);
}
/************************************************************************************/
void Foam::foamYade::initParticleForce(yadeParticle* particle) {

  vector nullV(0.0,0.0,0.0);
  particle->hydroForce = nullV;
  particle->hydroTorque = nullV;
}

/************************************************************************************/
void Foam::foamYade::calcHydroForce(){

  for (std::vector<yadeParticle*>::iterator pIter = localParticle.begin();  pIter != localParticle.end(); pIter++) {
    yadeParticle* particle = *pIter;
    initParticleForce(particle);
    if(isGaussianInterp){
      hydroDragForce(particle);
      buoyancyForce(particle);
    }
      else {
        particle -> hydroForce = vector(0,0,0);
        stokesDragForce(particle);
      }
   }
}

/************************************************************************************/
void Foam::foamYade::calcHydroTorque() {

    for (std::vector<yadeParticle*>::iterator pIter = localParticle.begin();  pIter != localParticle.end(); pIter++) {

      yadeParticle* particle = *pIter;

     if (isGaussianInterp){
       scalar s1 = 0.0; scalar s2 = 0.0; scalar s3 =0.0; 

     for (unsigned int i=0; i!= particle -> interpCellWeight.size(); ++i)
     {
       const int& cellid = particle -> interpCellWeight[i].first;
       const double& weight = particle -> interpCellWeight[i].second;
        s1 += ((vGrad[cellid].yz() - vGrad[cellid].zy())*weight);
        s2 += ((vGrad[cellid].zx() - vGrad[cellid].xz())*weight);
        s3 += ((vGrad[cellid].xy() - vGrad[cellid].yx())*weight); }

       const vector wfluid(s1,s2,s3);
       particle -> hydroTorque += M_PI*(pow(particle -> dia, 3))*(wfluid - particle-> rotationalVelocity)*nu*fluidDensity;

     } else {stokesDragTorque(particle);  }
  }
}
/************************************************************************************/
void Foam::foamYade::buoyancyForce(yadeParticle* particle) {

  vector bforce(0,0,0); double pv = 0.0;   const vector& g = gravity[0];

  for (unsigned int i=0; i != particle -> interpCellWeight.size(); ++i) {
    const double& wt =   particle -> interpCellWeight[i].second;
    pv += (particle->vol*wt);
  }
  bforce = (partDensity-fluidDensity)*pv*g;
  particle -> hydroForce += bforce;

  for (unsigned int i=0; i != particle -> interpCellWeight.size(); ++i ) {
    const double& weight = particle -> interpCellWeight[i].second;
    const int& id = particle -> interpCellWeight[i].first;
    const double& oo_cellVol = 1./(mesh.V()[id]*fluidDensity);
    uSource[id] = uSource[id] + (-bforce*weight*oo_cellVol);
  }
}

/************************************************************************************/

void Foam::foamYade::stokesDragForce(yadeParticle* particle) {

  autoPtr<interpolation<vector>> interpVel = interpolation<vector>::New("cell",U);  // cellPoint does not work in parallel, why?
  const vector& uFluid = interpVel->interpolate(particle->pos,particle->inCell);
  const double& coeff  = 3*M_PI*(particle->dia)*nu*fluidDensity;
  const double&  oo_cellVol = 1./(mesh.V()[particle-> inCell]*fluidDensity);
  particle->hydroForce = coeff*(uFluid-particle->linearVelocity);
  uSource[particle->inCell] +=  (-1*oo_cellVol*particle->hydroForce);
}

/************************************************************************************/

void Foam::foamYade::stokesDragTorque(yadeParticle* particle){

  autoPtr<interpolation<tensor>> interpGradU = interpolation<tensor>::New("cell",vGrad); //
  const tensor& uGradpt = interpGradU->interpolate(particle->pos, particle->inCell);
  scalar s1 = uGradpt.zy() - uGradpt.yz(); scalar s2 = uGradpt.zx()-uGradpt.xz(); scalar s3 = uGradpt.yx()-uGradpt.xy();
  vector wFluid(s1,s2,s3);
  particle->hydroTorque = M_PI *(pow(particle->dia,3))*(wFluid-particle->rotationalVelocity)*nu;
}


/*************************************************************************************/

void Foam::foamYade::setCellVolFraction(const std::vector<std::pair<label,double> >&  pVolContrib,
                                        const std::vector<std::pair<label, vector> > & uParticleContrib){

  if (pVolContrib.size()==0 || uParticleContrib.size()==0){return;}

  for (unsigned int i=0; i != pVolContrib.size(); ++i) {
    const label& id = pVolContrib[i].first;
    const double& pvolc = 1.0- (pVolContrib[i].second/mesh.V()[id]);
    alpha[id] = ((pvolc > 0) ? pvolc : 1e-08);
    const vector& upart = uParticleContrib[i].second;
    uParticle[id] = upart/mesh.V()[id];
  }

}
/************************************************************************************/

void Foam::foamYade::archimedesForce(yadeParticle* particle)
{
    // \vec F_a = (divT-divP)*pvol;
    vector divt(0,0,0); vector pg(0,0,0); double pv = 0.0;
    // get
    for (unsigned int i=0; i != particle -> interpCellWeight.size(); ++i) {

    const int& cellid = particle -> interpCellWeight[i].first;
    const double& weight = particle -> interpCellWeight[i].second;
    pv += (particle->vol*weight);
    divt = divt + (2.0*nu*divT[cellid]*weight);
    pg   =  pg + (-gradp[cellid]*weight);
  }

  //calculate
    const vector& f = pv*(divt+pg);
    particle -> hydroForce += f;

  // distribute
  for (unsigned int i=0; i != particle -> interpCellWeight.size(); ++i ) {
    const double& weight = particle -> interpCellWeight[i].second;
    const int& id = particle -> interpCellWeight[i].first;
    const double& oo_cellVol = 1./(mesh.V()[id]);
    uSource[id] = uSource[id] + (-f*weight*oo_cellVol);
  }

}
/************************************************************************************/

void Foam::foamYade::addedMassForce(yadeParticle* particle) {
    // const volVectorField& ddtU_f = fvc::ddt(U)+fvc::div(phi, U);
    // get
    vector ddtU(0,0,0); double pv = 0.0;
    for (unsigned int i=0; i != particle -> interpCellWeight.size(); ++i) {
      const double& weight = particle -> interpCellWeight[i].second;
      const label& cellid = particle -> interpCellWeight[i].first;
      pv += (particle->vol*weight);
      ddtU = ddtU + (ddtU_f[cellid]*weight);
  }
  //calculate
    const vector& f = pv*(ddtU - (particle->linearVelocity/deltaT))*partDensity;
    particle -> hydroForce += f;

  // distribute
  for (unsigned int i=0; i != particle -> interpCellWeight.size(); ++i ) {
    const double& weight = particle -> interpCellWeight[i].second;
    const int& id = particle -> interpCellWeight[i].first;
    const double& oo_cellVol = 1./(mesh.V()[id]);
    uSource[id] = uSource[id] + (-f*weight*oo_cellVol);
  }
}

/************************************************************************************/

void Foam::foamYade::hydroDragForce(yadeParticle* particle){

  if (particle -> interpCellWeight.size()==0) {return ; }

  vector uf(0,0,0);
  double alpha_p = 0;
  double small = 1e-09;
  double pv = 0.0;
  //get the velocities:
  //inteprolate fluid velocity to particle loc. and get the average particle vol frac

  for (unsigned int i=0; i != particle -> interpCellWeight.size(); ++i) {
    const int& cellid = particle->interpCellWeight[i].first;
    const double& wt = particle->interpCellWeight[i].second;
    uf = uf + (U[cellid]*wt);
    alpha_p = alpha_p + ((1-alpha[cellid])*wt);
    pv = pv + (particle->vol*wt);
  }

  alpha_p = alpha_p/(particle-> interpCellWeight.size());
  const double& alpha_f = 1-alpha_p;

  // calculate the force

  const vector& uRelVel = uf - particle->linearVelocity;
  const double& magUrelVel = mag(uRelVel);
  const scalar& Re = small+(magUrelVel*particle->dia)/nu;
  const double& pw = std::pow(Re, 0.687);
  const scalar& Cd = (24/Re)*(1+(0.15*pw));
  const double& coeff = (0.75*Cd*fluidDensity*magUrelVel)*(1/particle->dia)*(std::pow(alpha_f, -2.65));
  const vector& f = (pv)*coeff*uRelVel*alpha_f;

  particle-> hydroForce += f;

  //distribute the drag term "coeff" a.k.a 'K'. (negative coeff)

  for (unsigned int i=0; i != particle -> interpCellWeight.size(); ++i ) {
    const double& weight = particle -> interpCellWeight[i].second;
    const int& id = particle -> interpCellWeight[i].first;
    const scalar& value = (-coeff*weight*(1-alpha[id]));
    uSourceDrag[id] = uSourceDrag[id] + (value/fluidDensity);  //implicit source (scalar)
    uSource[id] = uSource[id] + ((value*uParticle[id])/fluidDensity); // explicit source term (with averaged particle velocity);
  }
}

/************************************************************************************/

void Foam::foamYade::combineForceToMaster(){

  double dummy;
  for (int i=0; i != numParticles; ++i) {
    std::vector<double> values(6, 0.0);
    haveIndex(i,values);
    for (int j=0; j<6; ++j){
      comm.procReduceSumDouble(values[j], dummy);
    }
  }
  resetLocalParticleList(localParticle);
}
/************************************************************************************/

void Foam::foamYade::setSourceZero() {

  forAll(uSource, cellI){
    uSource[cellI].x()=1e-15;
    uSource[cellI].y()=1e-15;
    uSource[cellI].z()=1e-15;
    if (isGaussianInterp){
        clearPvolContrib(pVolContrib, uParticleContrib);
        alpha[cellI] = 1.0;
        uSourceDrag[cellI] = 1e-15;
        uParticle[cellI].x() = 1e-15;
        uParticle[cellI].y() = 1e-15;
        uParticle[cellI].z() = 1e-15;

    }
  }
}
/************************************************************************************/
void Foam::foamYade::exchangeTimeStep() {
  if (comm.rank==1) {
    comm.sendOneDouble(yadeProc, deltaT);
    
  }
  double yadeDt;
  comm.cast_one_double(yadeProc, yadeDt);

}
/************************************************************************************/

void Foam::foamYade::haveIndex(const int& pIndex, std::vector<double>& values) {

  if (localParticle.size()==0){ return ;}

  for (unsigned int i=0; i != localParticle.size(); ++i) {
    if (localParticle[i]->indx==pIndex){
      values[0] = localParticle[i] -> hydroForce.x();
      values[1] = localParticle[i] -> hydroForce.y();
      values[2] = localParticle[i] -> hydroForce.z();
      values[3] = localParticle[i] -> hydroTorque.x();
      values[4] = localParticle[i] -> hydroTorque.y();
      values[5] = localParticle[i] -> hydroTorque.z();
    }
  }
}
/************************************************************************************/

void Foam::foamYade::sendHydroForceYade() {

  if(isGaussianInterp){
    combineForceToMaster();
  }else {
    sendHydroForcePoint();
  }
}

/************************************************************************************/
void Foam::foamYade::updateSources(std::vector<std::pair<label, double> >& pVolContrib ) {

  for (unsigned int i=0; i != pVolContrib.size(); ++i){

    const label& cellid  = pVolContrib[i].first;
    const double& pvolc = pVolContrib[i].second;
    const double& alpha_p = pvolc/mesh.V()[cellid];
    uSourceDrag[cellid] = alpha_p*uSourceDrag[cellid];
    uSource[cellid] = alpha_p*uSource[cellid];

  }
}

/************************************************************************************/


void Foam::foamYade::clearPvolContrib(std::vector<std::pair<label, double> >& pVolContrib,
  std::vector<std::pair<label,vector> >& uParticleContrib) {
  pVolContrib.clear();
  uParticleContrib.clear();
}


/************************************************************************************/

void Foam::foamYade::setParticleAction(double dt) {

  deltaT = dt;                      // get the time step
  locateAllParticle();              // find the particles in proc.
 if (isGaussianInterp){
    buildCelltoPartList(pVolContrib, uParticleContrib); // a mapping between the particles and it's associated cell ids -> sets the contribution of the particle volume and velocity
    setCellVolFraction(pVolContrib, uParticleContrib);  // calculate the fluid volfraction in the cells and set the particle velocity in the eulerian grid.
   }
  calcHydroForce();       // calculate the hydrodynamic forces
  calcHydroTorque();    // calculate the hydrodynamic torques.
  sendHydroForceYade();   // send info to yade
  exchangeTimeStep();     // exchange timesteps between yade and openfoam

}

/************************************************************************************/


void Foam::foamYade::recvTerminate(){
     comm.castTerminate();  
}
