/*  (c) 2019  Deepak Kunhappan  deepak.kn1990@gmail.com, deepak.kunhappan@3sr-grenoble.fr */ 

#ifndef FoamYadeMPI_H 
#define FoamYadeMPI_H 

#include <vector> 
#include <memory> 
#include <algorithm> 
#include <functional> 

#include "fvCFD.H" 
#include "meshTree.H" 
#include "interpolation.H" 
#include "interpolationCell.H" 
#include <mpi.h>

namespace Foam {
  
	class YadeParticle{
		public: 
			YadeParticle() {}; 
			double dia; 
			bool isFib = false; 
			double vol; 
			point pos;
			vector linearVelocity; 
			vector ori; 
			vector rotationalVelocity; 
			vector hydroForce; 
			vector hydroTorque; 
			std::vector<int> cellIds; 
			std::vector<std::pair<int, double> > interpCellWeight; 
			int inProc; 
			label inCell; 
			int indx; 
			void calcPartVol(){if (!isFib) vol = M_PI*std::pow(dia, 3.0)/6.0; }  
			bool sharedParticle = false; 
			std::vector<int> commRanks; 
			virtual ~YadeParticle(){}; 
			
	};
	
	class YadeProc{
		public : 
			YadeProc(){}; 
			int yRank;  							//  world rank of the yade Proc. 
			std::vector<int> numParticlesProc; 				// number of particles from intersections with yade proc
			bool inComm = false; 						//  flag for intersection with this proc. 
			std::vector<double> particleDataBuff;  				// buff to recv particle info -> pos, vel, dia, etc.
			std::vector<int> foundBuff; 	       				// buff to send grid search results, 
			std::vector<double> hydroForceBuff;    				// buff to send hydroforce and torque from each particle belonging to this proc.  
			std::vector<std::pair<int, double> > pVolContrib; 		// particle vol contribution from this yade proc. 
			std::vector<std::pair<int, vector> > uParticleContrib;  	// particle velocity contribution from this yade proc. 
			std::vector<std::shared_ptr<YadeParticle> > foundParticles;	// vector of found particles from this yadeProc. 
			int numParticles; 
			~YadeProc(){}; 
	}; 
  
	class FoamYade{
		protected: 
			/* tags for MPI messages */ 
			const int TAG_SZ_BUFF=1003; 
			const int TAG_GRID_BBOX = 1001; 
			const int TAG_YADE_DATA = 1002; 
			const int TAG_FORCE = 1005; 
			const int TAG_SEARCH_RES = 1004; 
			const int TAG_FLUID_DT = 1050; 
			const int TAG_YADE_DT = 1060; 
			const double small = 1e-09; 
	
		public: 
			std::vector<int> sendRanks; // used in serial to identify the force sending proc 
			std::vector<MPI_Request> reqVec; 
			int localRank, worldRank; 
			int localCommSize, worldCommSize; 
			int commSzDff; 
			bool rankSizeSet = false; 
			const fvMesh& mesh; 
			const volVectorField& U; 
			const volVectorField& gradP; 
			const volTensorField& vGrad; 
			const volVectorField& divT; 
			const volVectorField& ddtU; 
			const uniformDimensionedVectorField& g; 
			scalar rhoP; 
			scalar nu; 
			scalar rhoF; 
			bool isGaussianInterp; 
			volScalarField& uSourceDrag; 
			volScalarField& alpha; 
			volVectorField& uSource; 
			volVectorField& uParticle; 
			bool serialYade=false; 
			
			meshTree mshTree; 
			double yadeDT; 
			scalar deltaT; 
			scalar interpRange; 
			scalar interpRangeCu; 
			scalar sigmaInterp; 
			scalar sigmaPi; 
			std::vector<YadeProc> yadeProcs;
			std::vector<std::shared_ptr<YadeProc> > inCommProcs; 
			bool fibreCpl = false;  // flag for fiber coupling. 
			double hydorTimeScale; 
			double sigma; 
			
			FoamYade(const fvMesh& _mesh, 
				    const volVectorField& _U, 
				    const volVectorField& _gradP, 
				    const volTensorField& _vGrad, 
				    const volVectorField& _divT, 
				    const volVectorField& _ddtU, 
				    const uniformDimensionedVectorField& _g, 
				    volScalarField& _uSourceDrag, 
				    volScalarField& _alpha, 
				    volVectorField& _uSource, 
				    volVectorField& _uParticle, 
				    bool gaussianInterp) : 
				    mesh(_mesh), U(_U), gradP(_gradP), 
				    vGrad(_vGrad), divT(_divT), ddtU(_ddtU), 
				    g(_g), uSourceDrag(_uSourceDrag),  alpha(_alpha), 
				    uSource(_uSource), uParticle(_uParticle), mshTree(mesh.C())
				    {isGaussianInterp = gaussianInterp; getRankSize();  }
			void allocArrays(int, const std::shared_ptr<YadeProc>& ); 
			void getRankSize(); 
			void initFields(); 
			void setScalarProperties(scalar, scalar, scalar); 
			void sendMeshBbox(); 
			void recvYadeIntrs(); 
			void locateAllParticles(); 
			std::vector<int> locatePt(const vector& ); 
			void buildCellPartList(YadeProc* ); 
			void calcInterpWeightGaussian(std::vector<std::shared_ptr<YadeParticle> >& ); 
			void setCellVolFraction(YadeProc* ); 
			void updateSources(YadeProc* ); 
			void sendHydroForceYadeMPI(); 
			void exchangeDT();
			void sendHydroTimeScale(YadeProc* ); 
			void setParticleAction(double ); 
			void initParticleForce(YadeParticle* ); 
			void clearYadeProcs(); 
			void clearInCommProcs(); 
			//forces
			void calcHydroForce(YadeProc* ); 
			void calcHydroTorque(YadeProc* ); 
			// vol averaged 
			void hydroDragForce(YadeParticle*);
			void buoyancyForce(YadeParticle* ); 
			void addedMassForce(YadeParticle*); 
			void archimedesForce(YadeParticle* ); 
			void calcHydroTimeScale(); 
			// simple point force 
			void stokesDragForce(YadeParticle*); 
			void stokesDragTorque(YadeParticle*); 
			// clear sources 
			void setSourceZero(); 
			//for debug
			void printMsg(const std::string& ); 
			// terminate run 
			void finalizeRun(); 
			virtual ~FoamYade(){}; 
	}; 
  
}
#endif
