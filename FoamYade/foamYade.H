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
      scalar  vol;
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
      void calcPartVol(const double& dia){ vol =(M_PI*std::pow(dia,3.0))/6.0; }
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
       std::vector< std::pair <label, double> > pVolContrib;
       std::vector< std::pair <label, vector> > uParticleContrib;
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
       volScalarField& uSourceDrag;
       const volVectorField& ddtU_f;
       const uniformDimensionedVectorField& gravity;
       yadeComm comm;
       bool recvdParticleData;
       bool isGaussianInterp;
       meshTree mshTree;
       scalar partDensity;
       scalar nu;
       scalar fluidDensity;
       scalar deltaT;
       scalar interp_range;
       scalar interp_range_cu;
       scalar sigma_interp;
       scalar sigma_pi;

     protected:

     public:

       foamYade(const fvMesh& meshm,
                const volVectorField& u,
                      volVectorField& usource,
                      volVectorField& uparticle,
                      volScalarField& volfrac,
                const volVectorField& pgrad,
                const volTensorField& gradv,
                const volVectorField& divergencestress,
                      volScalarField& dragcoeff,
                const volVectorField& ddtuf,
                       uniformDimensionedVectorField& gravityf,
                       bool gaussianinterp) :
                      mesh(meshm), U(u), uSource(usource),
                      uParticle(uparticle), alpha(volfrac), gradp(pgrad),
                      vGrad(gradv), divT(divergencestress),
                      uSourceDrag(dragcoeff), ddtU_f(ddtuf),
                      gravity(gravityf), mshTree(mesh.C())
                      {allocArrays(); isGaussianInterp=gaussianinterp;  };
        void setScalarProperties(scalar, scalar, scalar);
        void locateAllParticle();
        bool locateParticle(yadeParticle* );
        void recvParticle();
        void sendHydroForcePoint();
        void setCellVolFraction(const std::vector<std::pair<label, double> >&,
                                const std::vector<std::pair<label, vector> >&);
        void initParticleForce(yadeParticle* );
        void calcInterpWeightGaussian(std::vector<yadeParticle*>&  );
        void calcInterpWeightDiracDelta(yadeParticle* );
        void archimedesForce(yadeParticle* );
        void hydroDragForce(yadeParticle* );
        void addedMassForce(yadeParticle* );
        void buoyancyForce(yadeParticle * );
        void calcHydroForce();
        void calcHydroTorque();
        void setSourceZero();
        void resetLocalParticleList(std::vector<yadeParticle*>& );
        void combineForceToMaster();
        void exchangeTimeStep();
        void sendHydroForceYade();
        void stokesDragForce(yadeParticle* );
        void stokesDragTorque(yadeParticle* );
        void haveIndex(const int&, std::vector<double>& );
        void buildCelltoPartList(std::vector<std::pair <label, double> >&,
                                std::vector<std::pair <label, vector> >& );
        void setParticleAction(double );
        void clearPvolContrib(std::vector<std::pair<label, double> >&,
                              std::vector<std::pair <label,vector> >& ); 

        void updateSources(std::vector<std::pair<label, double> >& ); 
        ~foamYade(){};

  };


}


#endif
