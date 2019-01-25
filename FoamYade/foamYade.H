#ifndef foamYade_H 
#define foamYade_H 

#include "fvCFD.H"
#include "yadeComm.H"
#include "interpolation.H"
#include "interpolationCell.H"
#include "interpolationCellPoint.H"
#include  <memory> 
#include  <vector> 

namespace Foam 
{ 

  class yadeParticle
  {
   
    public: 

      yadeParticle(){}; 
      scalar  dia; 
      scalar  reynoldsNum; 
      point   pos; 
      vector  linearVelocity; 
      vector  rotationalVelocity;
      vector  ori; 
      vector  hydroForce; 
      vector  hydroTorque;       
      int inProc;
      int indx; 
      label inCell;  
      ~yadeParticle(){}; 
  }; 


  class foamYade
  {
     private:
       
       int numParticles;  
       std::vector<double> hydroForce; 
       std::vector<double> particleData; 
       std::vector<int>  particleInProc;
       std::vector<yadeParticle*> localParticle;  
       bool haveParticle;
       const int yadeProc = 0; 
       yadeComm comm; 
       bool recvdParticleData;  



     public: 
       
       foamYade();
       void recvnumParticles(); 
       void locateAllParticle(fvMesh& mesh, scalar& nu_val,  volVectorField& U, volVectorField& uSource, volTensorField& gradU); 
       bool locateParticle(yadeParticle*,fvMesh& mesh); 
       void recvParticle();
       void sendHydroForce(); 
       void calcDragForce(fvMesh& mesh, yadeParticle*, scalar& nu ,volVectorField& U, volVectorField& uSource);
       void calcDragTorque(yadeParticle*, scalar& nu  , volTensorField& gradU); 
       void setSourceZero(volVectorField& uSource); 
       void resetlocalParticle();
       ~foamYade(){}; 

  };


}


#endif 
