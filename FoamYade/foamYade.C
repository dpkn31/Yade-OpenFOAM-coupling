//Deepak Kunhappan : deepak.kunhappan@3sr-grenoble.fr

#include "foamYade.H" 
#include <iostream> 


Foam::foamYade::foamYade()
 
{
  
  comm.cast_integer_data(yadeProc, numParticles);    // get the number of particles from YADE  
  particleData.reserve(numParticles*10); //pos, linearVel, rotVel, 'ori', dia;  
  particleInProc.reserve(numParticles);
  particleInProc.assign(numParticles,-1);
  hydroForce.assign(numParticles*6,1e-19); 
  haveParticle= false;
  recvdParticleData = false; 
    
}


void Foam::foamYade::resetlocalParticle()
{

  const std::size_t& vec_sz = localParticle.size(); 
  if (vec_sz > 0)
  {
    while (!localParticle.empty()) delete localParticle.back(), localParticle.pop_back(); 
  }
  localParticle.clear(); 

}



void Foam::foamYade::locateAllParticle(fvMesh& mesh, scalar& nu_val,volVectorField& U,
                                        volVectorField& uSource, 
                                        volTensorField& gradU) 
{
  

       comm.cast_double_array_data(yadeProc, particleData);      
       for (unsigned int i=0; i != (unsigned int) numParticles; ++i)
       {


         yadeParticle*  particle = new yadeParticle();
         //std::unique_ptr<yadeparticle-> particle (new yadeParticle); // TODO:use smartptr 
         particle->pos.x() = particleData[i*10];
         particle->pos.y()= particleData[i*10+1]; 
         particle->pos.z() = particleData[i*10+2]; 
         particle->inProc = -1; 

        
        if (locateParticle(particle, mesh)){
           
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

          calcDragForce(mesh, particle, nu_val, U, uSource); 
          calcDragTorque(particle,nu_val, gradU); 
          

        }

        comm.procReduce(particle->inProc,particleInProc[i]); 
        if (comm.rank==1 && particleInProc[i] < 0)
           std::cout << "particle->id = " << i  << " " << "not found " << std::endl;
        if (particle -> inProc != comm.rank)
          delete particle; 
      }
}


 
bool Foam::foamYade::locateParticle(yadeParticle* aParticle, fvMesh& mesh)    
{
  bool inside=false;  
  label value = mesh.findCell(aParticle->pos);
  if (value >= 0)
    inside =true; 
    aParticle->inCell = value;   
  return inside;            
 
}


void Foam::foamYade::sendHydroForce() 
{

  for (unsigned int i=0; i != (unsigned int) numParticles; ++i) 
  {
    if (comm.rank == particleInProc[i]) 
    {
      std::vector<yadeParticle*>::iterator pIter; 
        for (pIter=localParticle.begin(); pIter!= localParticle.end(); pIter++) 
        {
          if ((unsigned int) (*pIter)-> indx == i)
          {
            comm.sendOneDouble(yadeProc, (*pIter)->hydroForce.x());
            comm.sendOneDouble(yadeProc, (*pIter)->hydroForce.y());
            comm.sendOneDouble(yadeProc, (*pIter)->hydroForce.z());
            comm.sendOneDouble(yadeProc, (*pIter)->hydroTorque.x());
            comm.sendOneDouble(yadeProc, (*pIter)->hydroTorque.y());
            comm.sendOneDouble(yadeProc, (*pIter)->hydroTorque.z());
          }
        }
    }
  }
  resetlocalParticle(); 
}


void Foam::foamYade::setSourceZero(volVectorField& uSource)
{
  

  forAll(uSource, cellI)
  {
    uSource[cellI].x()=0.0; 
    uSource[cellI].y()=0.0; 
    uSource[cellI].z()=0.0; 
  }

} 


void Foam::foamYade::calcDragForce(fvMesh& mesh, yadeParticle* particle, scalar& nu_val,volVectorField& U, volVectorField& uSource) 
{

   //TODO: use Foam's drag laws. 
 
  autoPtr<interpolation<vector>> interpVel = interpolation<vector>::New("cell",U);  // cellPoint does not work in parallel, why? 
  const vector& uFluid = interpVel->interpolate(particle->pos,particle->inCell);
  particle->hydroForce = 3*M_PI*(particle->dia)*(uFluid-particle->linearVelocity)*nu_val; 
  uSource[particle->inCell] = -particle->hydroForce*(1/(mesh.V()[particle->inCell]));

   
}

void Foam::foamYade::calcDragTorque(yadeParticle* particle, scalar& nu_val, volTensorField& gradU) 

{

  //std::cout << "In calc torque " << std::endl;
  autoPtr<interpolation<tensor>> interpGradU = interpolation<tensor>::New("cell",gradU); //cellPoint doesn't work
  const tensor& uGradpt = interpGradU->interpolate(particle->pos, particle->inCell);
  scalar s1 = uGradpt.zy() - uGradpt.yz(); scalar s2 = uGradpt.zx()-uGradpt.xz(); scalar s3 = uGradpt.yx()-uGradpt.xy(); 
  vector wFluid(s1,s2,s3);  
  particle->hydroTorque = M_PI *(pow(particle->dia,3))*(wFluid-particle->rotationalVelocity)*nu_val;  
  
}
 

//void Foam::foamYade::calcRelaxationTime(); 
//separate class for fibers and pfacets? 




