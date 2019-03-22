//Deepak Kn , deepak.kunhappan@3sr-grenoble.fr, deepak.kn1990@gmail.com 

#include "foamYade.H" 
#include <iostream> 
#include <memory> 
#include <cmath>

void Foam::foamYade::allocArrays(){


  comm.cast_integer_data(yadeProc, numParticles); 
  hydroForce.reserve(numParticles*6);
  particleData.reserve(numParticles*10); // particleLinearVel, rotVel, ori, dia;  
  particleInProc.reserve(numParticles);
  particleInProc.assign(numParticles,-1);
  hydroForce.assign(numParticles*6,1e-19); 
  haveParticle= false;
  recvdParticleData = false; 
  numinCell.assign(mesh.V().size(), 0); 
  pVolContrib.assign(mesh.V().size(),0.0); 
  mshTree.build_tree();  
  interp_range = std::pow(mesh.V()[0], 0.3333)*3; // assuming uniform mesh, change later.
  sigma_interp = interp_range/2.3548200; // 2*deltaX/(2sqrt(2ln(2))) filter width half maximum; 
  sigma_pi = 1/(std::pow(2*M_PI*sigma_interp*sigma_interp, 1.5)); 
  
}


void Foam::foamYade::calcInterpWeightGaussian(yadeParticle* particle) { // gaussian distribution 

  if (! particle -> cellIds.size()) {return ; } 

  for (unsigned int i=0; i != particle -> cellIds.size(); ++i ) { 
    double distsq = 0; 
    const double& ds1 = mesh.C()[particle -> cellIds[i]].x() - particle -> pos.x(); 
    const double& ds2 = mesh.C()[particle -> cellIds[i]].y() - particle -> pos.y(); 
    const double& ds3 = mesh.C()[particle -> cellIds[i]].z() - particle -> pos.z(); 
    distsq = (ds1*ds1)+(ds2*ds2)+(ds3*ds3); 
    double weight = exp(-distsq/(2*std::pow(sigma_interp, 2)))*sigma_pi*mesh.V()[particle->cellIds[i]]; 

    particle -> interpCellWeight.push_back(std::make_pair(particle-> cellIds[i],weight)); 
    
  }

}


void Foam::foamYade::calcInterpWeightDiracDelta(yadeParticle* particle) { // gaussian distribution 

     std::cout << "Not implemented yet ! " << std::endl; 
//  if (! particle -> cellIds.size()) {return ; } 
//  for (unisgned int i=0; i != particle -> cellIds.size(); ++i ) { 
//    double distsq = 0; 
//    const double& ds1 = mesh.C()[particle -> cellIds[i]].x() - particle -> pos.x(); 
//    const double& ds2 = mesh.C()[particle -> cellIds[i]].y() - particle -> pos.y(); 
//    const double& ds3 = mesh.C()[particle -> cellIds[i]].z() - particle -> pos.z(); 
//    distsq = pow(ds1, 2)+pow(ds2, 2)+pow(ds3, 2);  
//    double weight = exp(-distsq/pow(sigma_interp, 2)) * sigma_p; 
//    particle -> interpCellWeight.push_back(std::make_pair(particle-> cellIds[i],weight)); 
//    
//  }
//
}




void Foam::foamYade::setScalarProperties(scalar nu_val, scalar pDensity, scalar fDensity) {
  nu = nu_val; 
  partDensity = pDensity; 
  fluidDensity = fDensity; 

}



void Foam::foamYade::resetLocalParticle()
{

  const std::size_t& vec_sz = localParticle.size(); 
  if (vec_sz > 0){
    while (!localParticle.empty()) delete localParticle.back(), localParticle.pop_back(); }
  localParticle.clear(); 

}



void Foam::foamYade::locateAllParticle() 
{
 
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
           particle->dia = particleData[i*10+9]; 
           setCellVolFraction(particle); 
           calcHydroForce(particle); 
           calcHydroTorque(particle); 
 
        }
        comm.procReduceMaxInt(particle-> inProc, particleInProc[i]);  
        if (comm.rank==1 && particleInProc[i] < 0)
           std::cout << "particle->id = " << i  << " " << " pos " << particle -> pos.x() << " " << particle -> pos.y() << " "<< particle -> pos.z() << " "<< "not found " << std::endl;
        if (particle -> inProc != comm.rank)
          delete particle; 
      }
}


 
bool Foam::foamYade::locateParticle(yadeParticle* aParticle)    
{

    bool value = false;
   if (isGaussianInterp ){ 
    aParticle ->cellIds = mshTree.nnearestCellsRange(aParticle->pos, interp_range, isGaussianInterp);
    if (aParticle->cellIds.size()){value = true; } 
   } else {  
     label cellid = mesh.findCell(aParticle -> pos); 
     if (cellid > -1 ) { aParticle->inCell= cellid; value = true;  }
   }
    return value; 
}



void Foam::foamYade::sendHydroForcePoint() 
{
 
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
  
  resetLocalParticle(); 
}


void Foam::foamYade::initParticleForce(yadeParticle* particle) {

  vector nullV(0.0,0.0,0.0); 
  particle->hydroForce = nullV; 
  particle->hydroTorque = nullV; 
}



void Foam::foamYade::calcHydroForce(yadeParticle* particle){

if(isGaussianInterp){

  initParticleForce(particle);
//   archimedesForce(particle); 
  hydroDragForce(particle); 
 // pressureForce(particle);  already included in archimedesForce; 
  // addedMassForce(particle); 
 //  buoyancyForce(particle); 
}
  else { 
    stokesDragForce(particle); 
  
  }

}


void Foam::foamYade::buoyancyForce(yadeParticle* particle) { 
 
  const scalar& pvol = M_PI*std::pow(particle -> dia, 3)*(1.0/6.0); 
  for (unsigned int i=0; i != particle -> interpCellWeight.size(); ++i) { 
  
    const int& cellid = particle -> interpCellWeight[i].first; 
    const double& weight = particle -> interpCellWeight[i].second; 
    const double& cellvol =  mesh.V()[cellid]; 
    const vector& vecgravity = gravity[cellid]*weight; 
    vector f = (partDensity - fluidDensity)*vecgravity*pvol; 
    uSource[cellid] += (-f/cellvol); 
    particle -> hydroForce += f; 
  
  }


}


void Foam::foamYade::stokesDragForce(yadeParticle* particle) {
  autoPtr<interpolation<vector>> interpVel = interpolation<vector>::New("cell",U);  // cellPoint does not work in parallel, why? 
  const vector& uFluid = interpVel->interpolate(particle->pos,particle->inCell);
  particle->hydroForce = 3*M_PI*(particle->dia)*(uFluid-particle->linearVelocity)*nu; 
  uSource[particle->inCell] = -particle->hydroForce*(1/(mesh.V()[particle->inCell]));
}

void Foam::foamYade::stokesDragTorque(yadeParticle* particle)  {

  autoPtr<interpolation<tensor>> interpGradU = interpolation<tensor>::New("cell",vGrad); //
  const tensor& uGradpt = interpGradU->interpolate(particle->pos, particle->inCell);
  scalar s1 = uGradpt.zy() - uGradpt.yz(); scalar s2 = uGradpt.zx()-uGradpt.xz(); scalar s3 = uGradpt.yx()-uGradpt.xy(); 
  vector wFluid(s1,s2,s3);  
  particle->hydroTorque = M_PI *(pow(particle->dia,3))*(wFluid-particle->rotationalVelocity)*nu;  

}



void Foam::foamYade::setCellVolFraction(yadeParticle* particle) //interpCellWeight vector<pair<id,weight>>
{

 if (isGaussianInterp){ 

    calcInterpWeightGaussian(particle);  
    scalar pvol = M_PI*pow(particle->dia,3)/6.0;
    for (unsigned int i=0; i != particle ->interpCellWeight.size(); ++i){
    const int& cellid = particle->interpCellWeight[i].first; 
    const double& weight = particle-> interpCellWeight[i].second; 
    const scalar& pf = (pvol* weight); 
    pVolContrib[cellid] += pf; 
    numinCell[cellid] +=1; 
    alpha[cellid] = 1-(pVolContrib[cellid]*numinCell[cellid])/mesh.V()[cellid]; 
  }
 
 } else { return ;}

}


void Foam::foamYade::archimedesForce(yadeParticle* particle)
{
    // \vec F_a = (divT-divP)*pvol; 
    const scalar& pvol =   2*nu*M_PI*pow(particle->dia,3)/6.0;
    for (unsigned int i=0; i != particle -> interpCellWeight.size(); ++i) {

    const int& cellid = particle -> interpCellWeight[i].first; 
    const double& weight = particle -> interpCellWeight[i].second; 
    const double& cellvol =  mesh.V()[cellid]; 

    const vector& divt  = 2.0*nu*divT[cellid]*weight;
    const vector& pg   =  -gradp[cellid]*weight; 
    vector f = pvol*(divt+pg); 
    uSource[cellid]+= -1*f/cellvol;   
    particle -> hydroForce += f; // fluid density 
  }

}



void Foam::foamYade::calcHydroTorque(yadeParticle* particle)
{

  if (isGaussianInterp){

  for (unsigned int i=0; i!= particle -> interpCellWeight.size(); ++i) 
  {
    const int& cellid = particle -> interpCellWeight[i].first; 
    const double& weight = particle -> interpCellWeight[i].second; 
    // const double& cellvol =  mesh.V()[cellid];  

    const scalar& s1 = (vGrad[cellid].yz() - vGrad[cellid].zy())*weight; 
    const scalar& s2 = (vGrad[cellid].zx() - vGrad[cellid].xz())*weight; 
    const scalar& s3 = (vGrad[cellid].xy() - vGrad[cellid].yx())*weight; 
    const vector wfluid(s1,s2,s3); 
    particle -> hydroTorque += M_PI*(pow(particle -> dia, 3))*(wfluid - particle-> rotationalVelocity)*nu*fluidDensity; 
  
  }

 } else {stokesDragTorque(particle);  }

}



void Foam::foamYade::addedMassForce(yadeParticle* particle) 
{
 // std::cout << "added mass not implemented yet" << std::endl;  
}


void Foam::foamYade::hydroDragForce(yadeParticle* particle){  

  // tp = rho_p(d**2)/18mu;                     
  const scalar& mu_ = nu*fluidDensity;
  const scalar& relaxTime = partDensity*std::pow(particle->dia, 2)/(18*mu_);
  for (unsigned int i=0; i != particle -> interpCellWeight.size(); ++i) {
  
    const int& cellid = particle -> interpCellWeight[i].first; 
    const double& weight = particle -> interpCellWeight[i].second; 
    const double& cellvol =  mesh.V()[cellid]; 
    const vector& ufluid = U[cellid]*weight;  

    const vector& rv = ufluid-particle->linearVelocity; 
    const scalar& Re = (particle->dia*mag(rv))/nu; 
    const vector& f  = 3*M_PI*mu_*rv*particle -> dia; 
    uSource[cellid]+= -1*f/cellvol; 
    vector& upart = uParticle[cellid]; 
    upart = upart + (weight*particle->linearVelocity); 
    particle -> hydroForce += f; 

  }

}


void Foam::foamYade::combineForceToMaster(){
  
  double dummy;
  for (int i=0; i != numParticles; ++i) { 
    std::vector<double> values(6, 1e-40);  
    haveIndex(i,values); 
    for (int j=0; j<6; ++j){
      comm.procReduceSumDouble(values[j], dummy); 
    }
  }
  resetLocalParticle(); 
}


void Foam::foamYade::setSourceZero() {
  
  forAll(uSource, cellI){
    uSource[cellI].x()=0.0; 
    uSource[cellI].y()=0.0; 
    uSource[cellI].z()=0.0; 
    alpha[cellI] = 1.0; 
    uParticle[cellI].x() = 0.0;
    uParticle[cellI].y() = 0.0; 
    uParticle[cellI].z() = 0.0; 
    std::fill(numinCell.begin(), numinCell.end(), 0); 
    std::fill(pVolContrib.begin(), pVolContrib.end(), 0.0); 
  } 
} 

void Foam::foamYade::exchangeTimeStep() { 
  if (comm.rank==1) { 
    comm.sendOneDouble(yadeProc, deltaT); 
  }
  double yadeDt =  1e-20;
  comm.cast_one_double(yadeProc, yadeDt); 

}


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

void Foam::foamYade::sendHydroForceYade() {

  if(isGaussianInterp){
    combineForceToMaster(); 
  }else {
    sendHydroForcePoint(); 
  }

}

void Foam::foamYade::setParticleAction(double dt) {

  deltaT = dt;
  locateAllParticle(); 
  sendHydroForceYade(); 
  exchangeTimeStep();  

}

