#ifndef FoamYadeMPI_H 
#define FoamYadeMPI_H 

#include "fvCFD.H" 
#include "meshTree.H" 
#include "interpolation.H" 
#include "interpolationCell.H" 
#include <vector> 
#include <memory> 
#include <algorithm> 
#include <functional> 

namespace Foam {
  
	class YadeParticle{
		public: 
			YadeParticle() {}; 
			double dia; 
			bool isFib; 
			double vol; 
			point pos;
			vector linearVelocity; 
			vector ori; 
			vector rotationalVelocity; 
			vector hydroForce; 
			vector hydroTroque; 
			std::vector<int> cellIds; 
			std::vector<std::pair<int, double> > interpCellWeight; 
			int inProc; 
			label inCell; 
			int indx; 
			void calcPartVol(){vol = M_PI*std::pow(dia, 3.0)/6.0; } 
			~YadeParticle(){}; 
			
	};
	
	class YadeProc{
		public : 
			YadeProc(){}; 
			int yRank;  							//  world rank of the yade Proc. 
			int numParticles; 						// number of particles from proc.
			bool inComm; 							//  flag for intersection with this proc. 
			std::vector<double> particleDataBuff;  				// buff to recv particle info -> pos, vel, dia, etc.
			std::vector<int> foundBuff; 	       				// buff to send grid search results, 
			std::vector<double> hydroForceBuff;    				// buff to send hydroforce and torque from each particle belonging to this proc.  
			std::vector<std::shared_ptr<YadeParticle> > particles;		// vector of 'found particles'
			std::vector<std::pair<int, double> > pVolContrib; 		// particle vol contribution from this yade proc. 
			std::vector<std::pair<int, vector> > uParticleContrib;  	// particle velocity contribution from this yade proc. 
			std::vector<std::shared_ptr<YadeParticle> > foundParticles;	// vector of found particles from this yadeProc. 
			~YadeProc(){}; 
	}; 
  
	class FoamYadeMPI{
		protected: 
			/* tags for MPI messages */ 
			const int TAG_SZ_BUFF=1003; 
			const int TAG_GRID_BBOX = 1001; 
			const int TAG_YADE_DATA = 1002; 
			const int TAG_FORCE = 1005; 
			const int TAG_SEARCH_RES = 1004; 
			const int TAG_FLUID_DT = 1050; 
			const int TAG_YADE_DT = 1060; 
	
		public: 
			int localRank, worldRank; 
			int localCommSize, worldCommSize; 
			int commSzDff; 
			bool rankSizeSet; 
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
			//yadeComm comm; 
			meshTree mshTree; 
			double yadeDT; 
			scalar deltaT; 
			scalar interpRange; 
			scalar interpRangeCu; 
			scalar sigmaInterp; 
			scalar sigmaPi; 
			std::vector<YadeProc> yadeProcs;
			std::vector<std::shared_ptr<YadeProc> > inCommProcs; 
			bool fibreCpl;  // flag for fiber coupling. 
			double hydorTimeScale; 
			double sigma
			
			FoamYadeMPI(const fvMesh& _mesh, 
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
				    scalar _rhoF, 
				    scalar _rhoP,
				    scalar _nu, 
				    bool gaussianInterp) : 
				    mesh(_mesh), U(_U), gradP(_gradP), 
				    vGrad(_vGrad), divT(_divT), uSourceDrag(_uSourceDrag), 
				    ddtU(_ddtU), g(_g),  alpha(_alpha), 
				    uSource(_uSource), uParticle(_uParticle), mshTree(mesh.C())
				    {isGaussianInterp = gaussianInterp; getRankSize(); setScalarProperties(_rhoP, _rhoF, _nu); }
			
			void getRankSize(); 
			void setScalarProperties(scalar, scalar, scalar); 
			void sendMeshBbox(); 
			void recvYadeIntrs(); 
			void locateAllParticles(); 
			void locatePt(const std::shared_ptr<YadeParticle>& ); 
			void buildCellPartList(std::shared_ptr<YadeProc>& ); 
			void calcInterpWeightGaussian(std::vector<std::shared_ptr<YadeParticle> >& ); 
			void setCellVolFraction(std::shared_ptr<YadeProc>& ); 
			void updateSources(std::shared_ptr<YadeProc>& ); 
			void sendHydroForceYadeMPI(); 
			void runCpl(); 
			void exchangeDT();
			void sendHydroTimeScale(); 
			void setParticleAction(const double& ); 
			void terminateRun(); 
			void initParticleForce(std::shared_ptr<YadeParticle>& ); 
			//forces
			void calcHydroForce(std::shared_ptr<YadeProc>& ); 
			void calcHydroTorque(std::shared_ptr<YadeParticle>& ); 
			// vol averaged 
			void hydroDragForce(std::shared_ptr<YadeParticle>& );
			void buoyancyForce(std::shared_ptr<YadeParticle>& ); 
			void addedMassForce(std::shared_ptr<YadeParticle>& ); 
			void archimedesForce(std::shared_ptr<YadeParticle>& ); 
			void calcHydroTimeScale(); 
			// simple point force 
			void stokesDragForce(std::shared_ptr<YadeParticle>& ); 
			void stokesDragTorque(std::shared_ptr<YadeParticle>& ); 
			~FoamYadeMPI(){}; 
	}; 
  
}
#endif