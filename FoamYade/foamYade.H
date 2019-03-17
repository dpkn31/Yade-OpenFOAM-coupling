// Deepak kunhappan, deepak.kunhappan@3sr-grenoble.fr / deepak.kn1990@gmail.com

#ifndef foamYade_H 
#define foamYade_H 

#include "fvCFD.H"
#include "yadeComm.H"
#include "meshTree.H"
#include "interpolation.H"
#include "interpolationCell.H"
#include "interpolationCellPoint.H"
#include  <vector> 
#include <utility> 
#include <algorithm> 
#include <functional> 

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
      std::vector<int> cellIds; 
      std::vector< std::pair<int, double> > interpCellWeight; 
      int inProc;
      int indx; 
      label inCell;
      long int nIter=0; 
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
       void allocArrays(); 
       const int yadeProc = 0;
       const fvMesh& mesh; 
       const volVectorField& U; 
       volVectorField& uSource;
       volVectorField& uParticle;  
       volScalarField& alpha; 
       const volVectorField& gradp;// fvc::grad(p) 
       const volTensorField& vGrad;// fvc::grad(U) 
       const volVectorField& divT; // fvc::laplacian((U)) 
       const uniformDimensionedVectorField& gravity; 
       yadeComm comm;
       bool recvdParticleData; 
       bool isGaussianIntep;  
       meshTree mshTree;
       scalar partDensity; 
       scalar nu;
       scalar fluidDensity; 
       scalar deltaT; 
       scalar interp_range; 
       scalar sigma_interp; 
       scalar sigma_pi; 

     protected:

     public: 
       
       foamYade(const fvMesh& meshm, 
                const volVectorField& u, 
                      volVectorField& usource, 
                      volVectorField& upart, 
                      volScalarField& volfrac,
                const volVectorField& pgrad, 
                const volTensorField& gradv,
                const volVectorField& divergencestress,  
                       uniformDimensionedVectorField&  gravityf, 
                       bool gaussianinterp) :
                      mesh(meshm), U(u), uSource(usource), 
                      uParticle(upart), alpha(volfrac), gradp(pgrad),
                      vGrad(gradv), divT(divergencestress), 
                      gravity(gravityf), mshTree(mesh.C()) 
                     {allocArrays(); isGaussianIntep=gaussianinterp;  }; 
        void setScalarProperties(scalar, scalar, scalar); 
        void locateAllParticle(); 
        bool locateParticle(yadeParticle* ); 
        void recvParticle();
        void sendHydroForcePoint(); 
        void setCellVolFraction(yadeParticle* );
        void initParticleForce(yadeParticle* );  
        void calcInterpWeightGaussian(yadeParticle* );
        void calcInterpWeightDiracDelta(yadeParticle* ); 
        void archimedesForce(yadeParticle* );  
        void hydroDragForce(yadeParticle* ); 
        void addedMassForce(yadeParticle* ); 
        void buoyancyForce(yadeParticle * ); 
        void calcHydroForce(yadeParticle* ); 
        void calcHydroTorque(yadeParticle* ); 
        void setSourceZero();  
        void resetLocalParticle(); 
        void combineForceToMaster(); 
        void exchangeTimeStep();  
        void sendHydroForceYade(); 
        void haveIndex(const int&, std::vector<double>& );
        void setParticleAction(double );  
        ~foamYade(){}; 

  };


}


#endif 
