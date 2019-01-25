#include "_utils.hpp"


#include <py/_utils.hpp>
#include <numpy/arrayobject.h>



static std::vector<int> bodyList; 
static std::vector<double> particleData; 
static std::vector<int> procList; 
static std::vector<double> hydroForce; 
static int rank, commSize;
const  int sendTag=500; 
static MPI_Status status;


//TODO: write in separate class, avoid static, try and catch if vecs not allocated.
//TODO: more abstraction to handle any type of particle.. 

void initMPI()
{

  //modify for YADE-MPI, 
  
  int threads; 
  MPI_Init_thread(0,0,MPI_THREAD_SINGLE,&threads); 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::cout  << "myrank = " << rank << std::endl;  
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);
  std::cout << "commSize = " << commSize << endl;  
  int color = 2; //Foam uses 1  
  MPI_Comm yadeComm; // dummy communicator; 
  MPI_Comm_split(MPI_COMM_WORLD,color,rank,&yadeComm);

}

void killMPI()
{
  MPI_Finalize(); 
}


void setnumParticles(int nparticles)

{
 bodyList.assign(nparticles, -1);
 particleData.assign(nparticles*10, 1e-19); 
 procList.assign(nparticles, -1);
 hydroForce.assign(nparticles*6, 1e-19);
 castNumpParticle(nparticles); 
}

void setBodylist(const py::list& abodylist)
{
 
  for (unsigned int i=0; i!= bodyList.capacity(); ++i)
  {
     bodyList[i] = py::extract<int> (abodylist[i]); 
  }

}

void castParticle() //TODO:: add OpenMP directives 
{
 shared_ptr<Scene> rb = Omega::instance().getScene(); 
 for (unsigned int i=0; i != bodyList.capacity() ; ++i)
 {
   const Body* b = (*rb -> bodies)[bodyList[i]].get(); 
   particleData[i*10+3] = b->state->vel[0]; 
   particleData[i*10+4] = b->state->vel[1]; 
   particleData[i*10+5] = b->state->vel[2]; 
   particleData[i*10+6] = b->state->angVel[0]; 
   particleData[i*10+7] = b->state->angVel[1]; 
   particleData[i*10+8] = b->state->angVel[2];
   shared_ptr<Sphere> s=YADE_PTR_DYN_CAST<Sphere>(b->shape);
   particleData[i*10+9] = s->radius; 
  // Vector3r pos;  
   if (rb -> isPeriodic){
     Vector3r pos = rb->cell->wrapPt(b->state->pos);
     particleData[i*10]   = pos[0];
     particleData[i*10+1] = pos[1];
     particleData[i*10+2] = pos[2];
   }

   else{
    
   particleData[i*10]   = b->state->pos[0]; 
   particleData[i*10+1] = b->state->pos[1]; 
   particleData[i*10+2] = b->state->pos[2]; 
   }
 
 }
   MPI_Bcast(&particleData.front(), particleData.capacity(), MPI_DOUBLE, rank, MPI_COMM_WORLD); 
}


void castNumpParticle(int& value) 
{
  
   MPI_Bcast(&value, 1, MPI_INT, rank, MPI_COMM_WORLD); 
 
} 


void updateProclist()
{
  for (unsigned int i=0; i != bodyList.capacity(); ++i) 
  {
     int dummy_val = -5; 
     MPI_Allreduce(&dummy_val,&procList[i],1,MPI_INT, MPI_MAX, MPI_COMM_WORLD);
     if (procList[i] < 0 ) 
       std::cout << "Particle not found in FOAM " << std::endl; 
   } 
}

void recvHydroForce() 
{
 
  for (unsigned int i=0; i!= procList.capacity(); ++i) 
 {
    int recvFrom = procList[i]; 
   for (unsigned int j=0; j != 6; ++j)
   {
     MPI_Recv(&hydroForce[6*i+j],1,MPI_DOUBLE,recvFrom,sendTag,MPI_COMM_WORLD,&status); 
   }
   
 }

}

void resetProcList()
{
   int capacity = procList.capacity(); 
   procList.assign(capacity, -1); 

}

void setHydroForce() 
{
   shared_ptr<Scene> rb = Omega::instance().getScene(); 
   for (unsigned int i=0; i!=bodyList.capacity(); ++i) {
    const Vector3r& fx=Vector3r(hydroForce[6*i], hydroForce[6*i+1], hydroForce[6*i+2]); 
    const Vector3r& tx=Vector3r(hydroForce[6*i+3], hydroForce[6*i+4], hydroForce[6*i+5]);
    //get body 
    //const Body* b1 = (*rb -> bodies)[boyList[i]].get(); 
    rb->forces.addForce(bodyList[i], fx); 
    rb->forces.addTorque(bodyList[i], tx);
  }
}

void runCoupling()
{

  castParticle(); 
  updateProclist(); 
  recvHydroForce();

}

py::tuple negPosExtremeIds(int axis, Real distFactor){
	py::tuple extrema=Shop::aabbExtrema();
	Real minCoord=py::extract<double>(extrema[0][axis])(),maxCoord=py::extract<double>(extrema[1][axis])();
	py::list minIds,maxIds;
	FOREACH(const shared_ptr<Body>& b, *Omega::instance().getScene()->bodies){
		shared_ptr<Sphere> s=YADE_PTR_DYN_CAST<Sphere>(b->shape); if(!s) continue;
		if(b->state->pos[axis]-s->radius*distFactor<=minCoord) minIds.append(b->getId());
		if(b->state->pos[axis]+s->radius*distFactor>=maxCoord) maxIds.append(b->getId());
	}
	return py::make_tuple(minIds,maxIds);
}

py::tuple coordsAndDisplacements(int axis,py::tuple Aabb){
	Vector3r bbMin(Vector3r::Zero()), bbMax(Vector3r::Zero()); bool useBB=py::len(Aabb)>0;
	if(useBB){bbMin=py::extract<Vector3r>(Aabb[0])();bbMax=py::extract<Vector3r>(Aabb[1])();}
	py::list retCoord,retDispl;
	FOREACH(const shared_ptr<Body>&b, *Omega::instance().getScene()->bodies){
		if(useBB && !Shop::isInBB(b->state->pos,bbMin,bbMax)) continue;
		retCoord.append(b->state->pos[axis]);
		retDispl.append(b->state->pos[axis]-b->state->refPos[axis]);
	}
	return py::make_tuple(retCoord,retDispl);
}

void setRefSe3(){
	Scene* scene=Omega::instance().getScene().get();
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		b->state->refPos=b->state->pos;
		b->state->refOri=b->state->ori;
	}
	if(scene->isPeriodic){
		scene->cell->refHSize=scene->cell->hSize;
	}
}

Real PWaveTimeStep(){return Shop::PWaveTimeStep();};
Real RayleighWaveTimeStep(){return Shop::RayleighWaveTimeStep();};

py::tuple interactionAnglesHistogram(int axis, int mask, size_t bins, py::tuple aabb, bool sphSph, Real minProjLen){
	if(axis<0||axis>2) throw invalid_argument("Axis must be from {0,1,2}=x,y,z.");
	Vector3r bbMin(Vector3r::Zero()), bbMax(Vector3r::Zero()); bool useBB=py::len(aabb)>0; if(useBB){bbMin=py::extract<Vector3r>(aabb[0])();bbMax=py::extract<Vector3r>(aabb[1])();}
	Real binStep=Mathr::PI/bins; int axis2=(axis+1)%3, axis3=(axis+2)%3;
	vector<Real> cummProj(bins,0.);
	shared_ptr<Scene> rb=Omega::instance().getScene();
	FOREACH(const shared_ptr<Interaction>& i, *rb->interactions){
		if(!i->isReal()) continue;
		const shared_ptr<Body>& b1=Body::byId(i->getId1(),rb), b2=Body::byId(i->getId2(),rb);
		if(!b1->maskOk(mask) || !b2->maskOk(mask)) continue;
		if(useBB && !Shop::isInBB(b1->state->pos,bbMin,bbMax) && !Shop::isInBB(b2->state->pos,bbMin,bbMax)) continue;
		if (sphSph && ( !dynamic_cast<Sphere*>(b1->shape.get()) || !dynamic_cast<Sphere*>(b2->shape.get()) ) ) continue;
		GenericSpheresContact* geom=dynamic_cast<GenericSpheresContact*>(i->geom.get());
		if(!geom) continue;
		Vector3r n(geom->normal); n[axis]=0.; Real nLen=n.norm();
		if(nLen<minProjLen) continue; // this interaction is (almost) exactly parallel to our axis; skip that one
		Real theta=acos(n[axis2]/nLen)*(n[axis3]>0?1:-1); if(theta<0) theta+=Mathr::PI;
		int binNo=theta/binStep;
		cummProj[binNo]+=nLen;
	}
	py::list val,binMid;
	for(size_t i=0; i<(size_t)bins; i++){ val.append(cummProj[i]); binMid.append(i*binStep);}
	return py::make_tuple(binMid,val);
}

py::tuple bodyNumInteractionsHistogram(py::tuple aabb){
	Vector3r bbMin(Vector3r::Zero()), bbMax(Vector3r::Zero()); bool useBB=py::len(aabb)>0; if(useBB){bbMin=py::extract<Vector3r>(aabb[0])();bbMax=py::extract<Vector3r>(aabb[1])();}
	const shared_ptr<Scene>& rb=Omega::instance().getScene();
	vector<int> bodyNumIntr; bodyNumIntr.resize(rb->bodies->size(),0);
	int maxIntr=0;
	FOREACH(const shared_ptr<Interaction>& i, *rb->interactions){
		if(!i->isReal()) continue;
		const Body::id_t id1=i->getId1(), id2=i->getId2(); const shared_ptr<Body>& b1=Body::byId(id1,rb), b2=Body::byId(id2,rb);
		if((useBB && Shop::isInBB(b1->state->pos,bbMin,bbMax)) || !useBB) {
			if (b1->isClumpMember()) bodyNumIntr[b1->clumpId]+=1; //count bodyNumIntr for the clump, not for the member 
			else bodyNumIntr[id1]+=1;
		}
		if((useBB && Shop::isInBB(b2->state->pos,bbMin,bbMax)) || !useBB) {
			if (b2->isClumpMember()) bodyNumIntr[b2->clumpId]+=1; //count bodyNumIntr for the clump, not for the member 
			else bodyNumIntr[id2]+=1;
		}
		maxIntr=max(max(maxIntr,bodyNumIntr[b1->getId()]),bodyNumIntr[b2->getId()]);
		if (b1->isClumpMember()) maxIntr=max(maxIntr,bodyNumIntr[b1->clumpId]);
		if (b2->isClumpMember()) maxIntr=max(maxIntr,bodyNumIntr[b2->clumpId]);
	}
	vector<int> bins; bins.resize(maxIntr+1,0);
	for(size_t id=0; id<bodyNumIntr.size(); id++){
		const shared_ptr<Body>& b=Body::byId(id,rb);
		if (b) {
			if(bodyNumIntr[id]>0) bins[bodyNumIntr[id]]+=1;
			// 0 is handled specially: add body to the 0 bin only if it is inside the bb requested (if applicable)
			// otherwise don't do anything, since it is outside the volume of interest
			else if(((useBB && Shop::isInBB(b->state->pos,bbMin,bbMax)) || !useBB) && !(b->isClumpMember())) bins[0]+=1;
		}
	}
	py::list count,num;
	for(size_t n=0; n<bins.size(); n++){
		if(bins[n]==0) continue;
		num.append(n); count.append(bins[n]);
	}
	return py::make_tuple(num,count);
}

Vector3r inscribedCircleCenter(const Vector3r& v0, const Vector3r& v1, const Vector3r& v2)
{
	return Shop::inscribedCircleCenter(v0,v1,v2);
}
py::dict getViscoelasticFromSpheresInteraction(Real tc, Real en, Real es)
{
	shared_ptr<ViscElMat> b = shared_ptr<ViscElMat>(new ViscElMat());
	Shop::getViscoelasticFromSpheresInteraction(tc,en,es,b);
	py::dict d;
	d["kn"]=b->kn;
	d["cn"]=b->cn;
	d["ks"]=b->ks;
	d["cs"]=b->cs;
	return d;
}
/* reset highlight of all bodies */
void highlightNone(){
	FOREACH(const shared_ptr<Body>& b, *Omega::instance().getScene()->bodies){
		if(!b->shape) continue;
		b->shape->highlight=false;
	}
}

/*!Sum moments acting on given bodies
 *
 * @param ids is the calculated bodies ids
 * @param axis is the direction of axis with respect to which the moment is calculated.
 * @param axisPt is a point on the axis.
 *
 * The computation is trivial: moment from force is is by definition r×F, where r
 * is position relative to axisPt; moment from moment is m; such moment per body is
 * projected onto axis.
 */
Real sumTorques(py::list ids, const Vector3r& axis, const Vector3r& axisPt){
	shared_ptr<Scene> rb=Omega::instance().getScene();
	rb->forces.sync();
	Real ret=0;
	size_t len=py::len(ids);
	for(size_t i=0; i<len; i++){
		const Body* b=(*rb->bodies)[py::extract<int>(ids[i])].get();
		const Vector3r& m=rb->forces.getTorque(b->getId());
		const Vector3r& f=rb->forces.getForce(b->getId());
		Vector3r r=b->state->pos-axisPt;
		ret+=axis.dot(m+r.cross(f));
	}
	return ret;
}
/* Sum forces acting on bodies within mask.
 *
 * @param ids list of ids
 * @param direction direction in which forces are summed
 *
 */
Real sumForces(py::list ids, const Vector3r& direction){
	shared_ptr<Scene> rb=Omega::instance().getScene();
	rb->forces.sync();
	Real ret=0;
	size_t len=py::len(ids);
	for(size_t i=0; i<len; i++){
		Body::id_t id=py::extract<int>(ids[i]);
		const Vector3r& f=rb->forces.getForce(id);
		ret+=direction.dot(f);
	}
	return ret;
}

/* Sum force acting on facets given by their ids in the sense of their respective normals.
   If axis is given, it will sum forces perpendicular to given axis only (not the the facet normals).
*/
Real sumFacetNormalForces(vector<Body::id_t> ids, int axis){
	shared_ptr<Scene> rb=Omega::instance().getScene(); rb->forces.sync();
	Real ret=0;
	FOREACH(const Body::id_t id, ids){
		Facet* f=YADE_CAST<Facet*>(Body::byId(id,rb)->shape.get());
		if(axis<0) ret+=rb->forces.getForce(id).dot(f->normal);
		else {
			Vector3r ff=rb->forces.getForce(id); ff[axis]=0;
			ret+=ff.dot(f->normal);
		}
	}
	return ret;
}

/* Set wire display of all/some/none bodies depending on the filter. */
void wireSome(string filter){
	enum{none,all,noSpheres,unknown};
	int mode=(filter=="none"?none:(filter=="all"?all:(filter=="noSpheres"?noSpheres:unknown)));
	if(mode==unknown) { LOG_WARN("Unknown wire filter `"<<filter<<"', using noSpheres instead."); mode=noSpheres; }
	FOREACH(const shared_ptr<Body>& b, *Omega::instance().getScene()->bodies){
		if(!b->shape) return;
		bool wire;
		switch(mode){
			case none: wire=false; break;
			case all: wire=true; break;
			case noSpheres: wire=!(bool)(YADE_PTR_DYN_CAST<Sphere>(b->shape)); break;
			default: throw logic_error("No such case possible");
		}
		b->shape->wire=wire;
	}
}
void wireAll(){wireSome("all");}
void wireNone(){wireSome("none");}
void wireNoSpheres(){wireSome("noSpheres");}


/* Tell us whether a point lies in polygon given by array of points.
 *  @param xy is the point that is being tested
 *  @param vertices is Numeric.array (or list or tuple) of vertices of the polygon.
 *         Every row of the array is x and y coordinate, numer of rows is >= 3 (triangle).
 *
 * Copying the algorithm from http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
 * is gratefully acknowledged:
 *
 * License to Use:
 * Copyright (c) 1970-2003, Wm. Randolph Franklin
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *   1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
 *   2. Redistributions in binary form must reproduce the above copyright notice in the documentation and/or other materials provided with the distribution.
 *   3. The name of W. Randolph Franklin may not be used to endorse or promote products derived from this Software without specific prior written permission. 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * http://numpy.scipy.org/numpydoc/numpy-13.html told me how to use Numeric.array from c
 */
bool pointInsidePolygon(py::tuple xy, py::object vertices){
	Real testx=py::extract<double>(xy[0])(),testy=py::extract<double>(xy[1])();
	char** vertData; int rows, cols; PyArrayObject* vert=(PyArrayObject*)vertices.ptr();
	int result=PyArray_As2D((PyObject**)&vert /* is replaced */ ,&vertData,&rows,&cols,PyArray_DOUBLE);
	if(result!=0) throw invalid_argument("Unable to cast vertices to 2d array");
	if(cols!=2 || rows<3) throw invalid_argument("Vertices must have 2 columns (x and y) and at least 3 rows.");
	int i /*current node*/, j/*previous node*/; bool inside=false;
	for(i=0,j=rows-1; i<rows; j=i++){
		double vx_i=*(double*)(vert->data+i*vert->strides[0]), vy_i=*(double*)(vert->data+i*vert->strides[0]+vert->strides[1]), vx_j=*(double*)(vert->data+j*vert->strides[0]), vy_j=*(double*)(vert->data+j*vert->strides[0]+vert->strides[1]);
		if (((vy_i>testy)!=(vy_j>testy)) && (testx < (vx_j-vx_i) * (testy-vy_i) / (vy_j-vy_i) + vx_i) ) inside=!inside;
	}
	Py_DECREF(vert);
	return inside;
}

/* Compute area of convex hull when when taking (swept) spheres crossing the plane at coord, perpendicular to axis.

	All spheres that touch the plane are projected as hexagons on their circumference to the plane.
	Convex hull from this cloud is computed.
	The area of the hull is returned.

*/
Real approxSectionArea(Real coord, int axis){
	std::list<Vector2r> cloud;
	if(axis<0 || axis>2) throw invalid_argument("Axis must be ∈ {0,1,2}");
	const int ax1=(axis+1)%3, ax2=(axis+2)%3;
	const Real sqrt3=sqrt(3);
	Vector2r mm,mx; int i=0;
	FOREACH(const shared_ptr<Body>& b, *Omega::instance().getScene()->bodies){
		Sphere* s=dynamic_cast<Sphere*>(b->shape.get());
		if(!s) continue;
		const Vector3r& pos(b->state->pos); const Real r(s->radius);
		if((pos[axis]>coord && (pos[axis]-r)>coord) || (pos[axis]<coord && (pos[axis]+r)<coord)) continue;
		Vector2r c(pos[ax1],pos[ax2]);
		cloud.push_back(c+Vector2r(r,0.)); cloud.push_back(c+Vector2r(-r,0.));
		cloud.push_back(c+Vector2r( r/2., sqrt3*r)); cloud.push_back(c+Vector2r( r/2.,-sqrt3*r));
		cloud.push_back(c+Vector2r(-r/2., sqrt3*r)); cloud.push_back(c+Vector2r(-r/2.,-sqrt3*r));
		if(i++==0){ mm=c, mx=c;}
		mm=Vector2r(min(c[0]-r,mm[0]),min(c[1]-r,mm[1]));
		mx=Vector2r(max(c[0]+r,mx[0]),max(c[1]+r,mx[1]));
	}
	if(cloud.size()==0) return 0;
	ConvexHull2d ch2d(cloud);
	vector<Vector2r> hull=ch2d();
	return simplePolygonArea2d(hull);
}
/* Find all interactions deriving from NormShearPhys that cross plane given by a point and normal
	(the normal may not be normalized in this case, though) and sum forces (both normal and shear) on them.
	
	Returns a 3-tuple with the components along global x,y,z axes, which can be viewed as "action from lower part, towards
	upper part" (lower and upper parts with respect to the plane's normal).

	(This could be easily extended to return sum of only normal forces or only of shear forces.)
*/
Vector3r forcesOnPlane(const Vector3r& planePt, const Vector3r&  normal){
	Vector3r ret(Vector3r::Zero());
	Scene* scene=Omega::instance().getScene().get();
	FOREACH(const shared_ptr<Interaction>&I, *scene->interactions){
		if(!I->isReal()) continue;
		NormShearPhys* nsi=dynamic_cast<NormShearPhys*>(I->phys.get());
		if(!nsi) continue;
		Vector3r pos1,pos2;
		pos1=Body::byId(I->getId1(),scene)->state->pos; pos2=Body::byId(I->getId2(),scene)->state->pos;
		Real dot1=(pos1-planePt).dot(normal), dot2=(pos2-planePt).dot(normal);
		if(dot1*dot2>0) continue; // both (centers of) bodies are on the same side of the plane=> this interaction has to be disregarded
		// if pt1 is on the negative plane side, d3dg->normal.Dot(normal)>0, the force is well oriented;
		// otherwise, reverse its contribution. So that we return finally
		// Sum [ ( normal(plane) dot normal(interaction= from 1 to 2) ) "nsi->force" ]
		ret+=(dot1<0.?1.:-1.)*(nsi->normalForce+nsi->shearForce);
	}
	return ret;
}

/* Less general than forcesOnPlane, computes force on plane perpendicular to axis, passing through coordinate coord. */
Vector3r forcesOnCoordPlane(Real coord, int axis){
	Vector3r planePt(Vector3r::Zero()); planePt[axis]=coord;
	Vector3r normal(Vector3r::Zero()); normal[axis]=1;
	return forcesOnPlane(planePt,normal);
}


py::tuple spiralProject(const Vector3r& pt, Real dH_dTheta, int axis, Real periodStart, Real theta0){
	Real r,h,theta;
	boost::tie(r,h,theta)=Shop::spiralProject(pt,dH_dTheta,axis,periodStart,theta0);
	return py::make_tuple(py::make_tuple(r,h),theta);
}

shared_ptr<Interaction> Shop__createExplicitInteraction(Body::id_t id1, Body::id_t id2){ return Shop::createExplicitInteraction(id1,id2,/*force*/true);}

Real Shop__unbalancedForce(bool useMaxForce /*false by default*/){return Shop::unbalancedForce(useMaxForce);}
py::tuple Shop__totalForceInVolume(){Real stiff; Vector3r ret=Shop::totalForceInVolume(stiff); return py::make_tuple(ret,stiff); }
Real Shop__getSpheresVolume(int mask){ return Shop::getSpheresVolume(Omega::instance().getScene(), mask);}
Real Shop__getSpheresMass(int mask){ return Shop::getSpheresMass(Omega::instance().getScene(), mask);}
py::object Shop__kineticEnergy(bool findMaxId){
	if(!findMaxId) return py::object(Shop::kineticEnergy());
	Body::id_t maxId;
	Real E=Shop::kineticEnergy(NULL,&maxId);
	return py::make_tuple(E,maxId);
}

Real maxOverlapRatio(){
	Scene* scene=Omega::instance().getScene().get();
	Real ret=-1;
	FOREACH(const shared_ptr<Interaction> I, *scene->interactions){
		if(!I->isReal()) continue;
		Sphere *s1(dynamic_cast<Sphere*>(Body::byId(I->getId1(),scene)->shape.get())), *s2(dynamic_cast<Sphere*>(Body::byId(I->getId2(),scene)->shape.get()));
		if((!s1) || (!s2)) continue;
		ScGeom* geom=dynamic_cast<ScGeom*>(I->geom.get());
		if(!geom) continue;
		Real rEq=2*s1->radius*s2->radius/(s1->radius+s2->radius);
		ret=max(ret,geom->penetrationDepth/rEq);
	}
	return ret;
}

Real Shop__getPorosity(Real volume){ return Shop::getPorosity(Omega::instance().getScene(),volume); }
Real Shop__getVoxelPorosity(int resolution, Vector3r start,Vector3r end){ return Shop::getVoxelPorosity(Omega::instance().getScene(),resolution,start,end); }

//Matrix3r Shop__stressTensorOfPeriodicCell(bool smallStrains=false){return Shop::stressTensorOfPeriodicCell(smallStrains);}
py::tuple Shop__fabricTensor(Real cutoff, bool splitTensor, Real thresholdForce){return Shop::fabricTensor(cutoff, splitTensor,thresholdForce);}
py::tuple Shop__normalShearStressTensors(bool compressionPositive, bool splitNormalTensor, Real thresholdForce){return Shop::normalShearStressTensors(compressionPositive,splitNormalTensor,thresholdForce);}

py::list Shop__getStressLWForEachBody(){return Shop::getStressLWForEachBody();}

py::list Shop__getBodyIdsContacts(Body::id_t bodyID){return Shop::getBodyIdsContacts(bodyID);}

Real shiftBodies(py::list ids, const Vector3r& shift){
	shared_ptr<Scene> rb=Omega::instance().getScene();
	size_t len=py::len(ids);
	for(size_t i=0; i<len; i++){
		const Body* b=(*rb->bodies)[py::extract<int>(ids[i])].get();
		if(!b) continue;
		b->state->pos+=shift;
	}
	return 1;
}

void Shop__calm(int mask){ return Shop::calm(Omega::instance().getScene(), mask);}

void setNewVerticesOfFacet(const shared_ptr<Body>& b, const Vector3r& v1, const Vector3r& v2, const Vector3r& v3) {
	Vector3r center = inscribedCircleCenter(v1,v2,v3);
	Facet* facet = YADE_CAST<Facet*>(b->shape.get());
	facet->vertices[0] = v1 - center;
	facet->vertices[1] = v2 - center;
	facet->vertices[2] = v3 - center;
	b->state->pos = center;
}

py::list intrsOfEachBody() {
	py::list ret, temp;
	shared_ptr<Scene> rb=Omega::instance().getScene();
	size_t n = rb->bodies->size();
	// create list of size len(O.bodies) full of zeros
	for (size_t i=0; i<n; i++) {
		ret.append(py::list());
	}
	// loop over all interactions and fill the list ret
	FOREACH(const shared_ptr<Interaction>& i, *rb->interactions) {
		if (!i->isReal()) { continue; }
		temp = py::extract<py::list>(ret[i->getId1()]);
		temp.append( i );
		temp = py::extract<py::list>(ret[i->getId2()]);
		temp.append( i );
	}
	return ret;
}

py::list numIntrsOfEachBody() {
	py::list ret;
	shared_ptr<Scene> rb=Omega::instance().getScene();
	size_t n = rb->bodies->size();
	// create list of size len(O.bodies) full of zeros
	for (size_t i=0; i<n; i++) {
		ret.append(0);
	}
	// loop over all interactions and fill the list ret
	FOREACH(const shared_ptr<Interaction>& i, *rb->interactions) {
		if (!i->isReal()) continue;
		ret[i->getId1()] += 1;
		ret[i->getId2()] += 1;
	}
	return ret;
}

Real Shop__getSpheresVolume2D(int mask=-1){ return Shop::getSpheresVolume2D(Omega::instance().getScene(), mask);}
Real Shop__getVoidRatio2D(Real zlen=1){ return Shop::getVoidRatio2D(Omega::instance().getScene(),zlen);}
py::tuple Shop__getStressAndTangent(Real volume=0, bool symmetry=true){return Shop::getStressAndTangent(volume,symmetry);}

void setBodyPosition(int id, Vector3r newPos, string axis){
	shared_ptr<Scene> rb=Omega::instance().getScene();
	const Body* b=(*rb->bodies)[id].get();
	for (char c : axis){
		if(c=='x') {b->state->pos[0] = newPos[0];continue;}
		if(c=='y') {b->state->pos[1] = newPos[1];continue;}
		if(c=='z') {b->state->pos[2] = newPos[2];continue;}
	}
}

void setBodyVelocity(int id, Vector3r newVel, string axis){
	shared_ptr<Scene> rb=Omega::instance().getScene();
	const Body* b=(*rb->bodies)[id].get();
	for (char c : axis){
		if(c=='x') {b->state->vel[0] = newVel[0];continue;}
		if(c=='y') {b->state->vel[1] = newVel[1];continue;}
		if(c=='z') {b->state->vel[2] = newVel[2];continue;}
	}
}

void setBodyOrientation(int id, Quaternionr newOri){
	shared_ptr<Scene> rb=Omega::instance().getScene();
	const Body* b=(*rb->bodies)[id].get();
	b->state->ori=newOri;
}

void setBodyAngularVelocity(int id, Vector3r newAngVel){
	shared_ptr<Scene> rb=Omega::instance().getScene();
	const Body* b=(*rb->bodies)[id].get();
	b->state->angVel=newAngVel;
}

void setBodyColor(int id, Vector3r newColor){
  shared_ptr<Scene> rb=Omega::instance().getScene();
	const Body* b=(*rb->bodies)[id].get();
	b->shape->color=newColor;
}

BOOST_PYTHON_MODULE(_utils){
	YADE_SET_DOCSTRING_OPTS;

  py::def("initMPI", initMPI, "initialize communicator for coupling with OpenFOAM.");
  py::def("killMPI", killMPI, "destory mpi environment"); 
  py::def("setHydroForce", setHydroForce, "set hydro force and torque on bodies");
  py::def("setnumParticles", setnumParticles, "number if particles ,");  
//  py::def("castParticle", castParticle, "MPI Broadcast of particle information to OpenFOAM procs."); 
//  py::def("updateProclist",updateProclist, "ranks which has particle in their local mesh and id to recieve from."); 
//  py::def("resetProcList",resetProcList, "re-initialize the procList array."); 
//  py::def("recvHydroForce", recvHydroForce, "receive the hydro force and torque from OpenFOAM"); 
  py::def("setBodylist", setBodylist, (py::arg("bodyList")),  "body ids particiipating in fluid interaction");
  py::def("runCoupling", runCoupling, "function for data coupling");  
	py::def("PWaveTimeStep",PWaveTimeStep,"Get timestep accoring to the velocity of P-Wave propagation; computed from sphere radii, rigidities and masses.");
	py::def("RayleighWaveTimeStep",RayleighWaveTimeStep,"Determination of time step according to Rayleigh wave speed of force propagation.");
	py::def("getSpheresVolume",Shop__getSpheresVolume,(py::arg("mask")=-1),"Compute the total volume of spheres in the simulation, mask parameter is considered");
	py::def("getSpheresMass",Shop__getSpheresMass,(py::arg("mask")=-1),"Compute the total mass of spheres in the simulation, mask parameter is considered");
	py::def("porosity",Shop__getPorosity,(py::arg("volume")=-1),"Compute packing porosity $\\frac{V-V_s}{V}$ where $V$ is overall volume and $V_s$ is volume of spheres.\n\n:param float volume: overall volume $V$. For periodic simulations, current volume of the :yref:`Cell` is used. For aperiodic simulations, the value deduced from utils.aabbDim() is used. For compatibility reasons, positive values passed by the user are also accepted in this case.\n");
	py::def("voxelPorosity",Shop__getVoxelPorosity,(py::arg("resolution")=200,py::arg("start")=Vector3r(0,0,0),py::arg("end")=Vector3r(0,0,0)),"Compute packing porosity $\\frac{V-V_v}{V}$ where $V$ is a specified volume (from start to end) and $V_v$ is volume of voxels that fall inside any sphere. The calculation method is to divide whole volume into a dense grid of voxels (at given resolution), and count the voxels that fall inside any of the spheres. This method allows one to calculate porosity in any given sub-volume of a whole sample. It is properly excluding part of a sphere that does not fall inside a specified volume.\n\n:param int resolution: voxel grid resolution, values bigger than resolution=1600 require a 64 bit operating system, because more than 4GB of RAM is used, a resolution=800 will use 500MB of RAM.\n:param Vector3 start: start corner of the volume.\n:param Vector3 end: end corner of the volume.\n");
	py::def("aabbExtrema",Shop::aabbExtrema,(py::arg("cutoff")=0.0,py::arg("centers")=false),"Return coordinates of box enclosing all spherical bodies\n\n:param bool centers: do not take sphere radii in account, only their centroids\n:param float∈〈0…1〉 cutoff: relative dimension by which the box will be cut away at its boundaries.\n\n\n:return: [lower corner, upper corner] as [Vector3,Vector3]\n\n");
	py::def("ptInAABB",Shop::isInBB,"Return True/False whether the point p is within box given by its min and max corners");
	py::def("negPosExtremeIds",negPosExtremeIds,(py::arg("axis"),py::arg("distFactor")),"Return list of ids for spheres (only) that are on extremal ends of the specimen along given axis; distFactor multiplies their radius so that sphere that do not touch the boundary coordinate can also be returned.");
	py::def("approxSectionArea",approxSectionArea,"Compute area of convex hull when when taking (swept) spheres crossing the plane at coord, perpendicular to axis.");
	py::def("coordsAndDisplacements",coordsAndDisplacements,(py::arg("axis"),py::arg("Aabb")=py::tuple()),"Return tuple of 2 same-length lists for coordinates and displacements (coordinate minus reference coordinate) along given axis (1st arg); if the Aabb=((x_min,y_min,z_min),(x_max,y_max,z_max)) box is given, only bodies within this box will be considered.");
	py::def("setRefSe3",setRefSe3,"Set reference :yref:`positions<State::refPos>` and :yref:`orientations<State::refOri>` of all :yref:`bodies<Body>` equal to their current :yref:`positions<State::pos>` and :yref:`orientations<State::ori>`.");
	py::def("interactionAnglesHistogram",interactionAnglesHistogram,(py::arg("axis"),py::arg("mask")=0,py::arg("bins")=20,py::arg("aabb")=py::tuple(),py::arg("sphSph")=0,py::arg("minProjLen")=1e-6));
	py::def("bodyNumInteractionsHistogram",bodyNumInteractionsHistogram,(py::arg("aabb")));
	py::def("inscribedCircleCenter",inscribedCircleCenter,(py::arg("v1"),py::arg("v2"),py::arg("v3")),"Return center of inscribed circle for triangle given by its vertices *v1*, *v2*, *v3*.");
	py::def("unbalancedForce",&Shop__unbalancedForce,(py::args("useMaxForce")=false),"Compute the ratio of mean (or maximum, if *useMaxForce*) summary force on bodies and mean force magnitude on interactions. For perfectly static equilibrium, summary force on all bodies is zero (since forces from interactions cancel out and induce no acceleration of particles); this ratio will tend to zero as simulation stabilizes, though zero is never reached because of finite precision computation. Sufficiently small value can be e.g. 1e-2 or smaller, depending on how much equilibrium it should be.");
	py::def("kineticEnergy",Shop__kineticEnergy,(py::args("findMaxId")=false),"Compute overall kinetic energy of the simulation as\n\n.. math:: \\sum\\frac{1}{2}\\left(m_i\\vec{v}_i^2+\\vec{\\omega}(\\mat{I}\\vec{\\omega}^T)\\right).\n\nFor :yref:`aspherical<Body.aspherical>` bodies, the inertia tensor $\\mat{I}$ is transformed to global frame, before multiplied by $\\vec{\\omega}$, therefore the value should be accurate.\n");
	py::def("sumForces",sumForces,(py::arg("ids"),py::arg("direction")),"Return summary force on bodies with given *ids*, projected on the *direction* vector.");
	py::def("sumTorques",sumTorques,(py::arg("ids"),py::arg("axis"),py::arg("axisPt")),"Sum forces and torques on bodies given in *ids* with respect to axis specified by a point *axisPt* and its direction *axis*.");
	py::def("sumFacetNormalForces",sumFacetNormalForces,(py::arg("ids"),py::arg("axis")=-1),"Sum force magnitudes on given bodies (must have :yref:`shape<Body.shape>` of the :yref:`Facet` type), considering only part of forces perpendicular to each :yref:`facet's<Facet>` face; if *axis* has positive value, then the specified axis (0=x, 1=y, 2=z) will be used instead of facet's normals.");
	py::def("forcesOnPlane",forcesOnPlane,(py::arg("planePt"),py::arg("normal")),"Find all interactions deriving from :yref:`NormShearPhys` that cross given plane and sum forces (both normal and shear) on them.\n\n:param Vector3 planePt: a point on the plane\n:param Vector3 normal: plane normal (will be normalized).\n");
	py::def("forcesOnCoordPlane",forcesOnCoordPlane);
	py::def("totalForceInVolume",Shop__totalForceInVolume,"Return summed forces on all interactions and average isotropic stiffness, as tuple (Vector3,float)");
	py::def("createInteraction",Shop__createExplicitInteraction,(py::arg("id1"),py::arg("id2")),"Create interaction between given bodies by hand.\n\nCurrent engines are searched for :yref:`IGeomDispatcher` and :yref:`IPhysDispatcher` (might be both hidden in :yref:`InteractionLoop`). Geometry is created using ``force`` parameter of the :yref:`geometry dispatcher<IGeomDispatcher>`, wherefore the interaction will exist even if bodies do not spatially overlap and the functor would return ``false`` under normal circumstances. \n\n.. warning:: This function will very likely behave incorrectly for periodic simulations (though it could be extended it to handle it farily easily).");
	py::def("spiralProject",spiralProject,(py::arg("pt"),py::arg("dH_dTheta"),py::arg("axis")=2,py::arg("periodStart")=std::numeric_limits<Real>::quiet_NaN(),py::arg("theta0")=0));
	py::def("pointInsidePolygon",pointInsidePolygon);
	py::def("scalarOnColorScale",Shop::scalarOnColorScale);
	py::def("highlightNone",highlightNone,"Reset :yref:`highlight<Shape::highlight>` on all bodies.");
	py::def("wireAll",wireAll,"Set :yref:`Shape::wire` on all bodies to True, rendering them with wireframe only.");
	py::def("wireNone",wireNone,"Set :yref:`Shape::wire` on all bodies to False, rendering them as solids.");
	py::def("wireNoSpheres",wireNoSpheres,"Set :yref:`Shape::wire` to True on non-spherical bodies (:yref:`Facets<Facet>`, :yref:`Walls<Wall>`).");
	py::def("flipCell",&Shop::flipCell,(py::arg("flip")=Matrix3r(Matrix3r::Zero())),"Flip periodic cell so that angles between $R^3$ axes and transformed axes are as small as possible. This function relies on the fact that periodic cell defines by repetition or its corners regular grid of points in $R^3$; however, all cells generating identical grid are equivalent and can be flipped one over another. This necessiatates adjustment of :yref:`Interaction.cellDist` for interactions that cross boundary and didn't before (or vice versa), and re-initialization of collider. The *flip* argument can be used to specify desired flip: integers, each column for one axis; if zero matrix, best fit (minimizing the angles) is computed automatically.\n\nIn c++, this function is accessible as ``Shop::flipCell``.");
	py::def("getViscoelasticFromSpheresInteraction",getViscoelasticFromSpheresInteraction,(py::arg("tc"),py::arg("en"),py::arg("es")),"Attention! The function is deprecated! Compute viscoelastic interaction parameters from analytical solution of a pair spheres collision problem:\n\n.. math:: k_n=\\frac{m}{t_c^2}\\left(\\pi^2+(\\ln e_n)^2\\right) \\\\ c_n=-\\frac{2m}{t_c}\\ln e_n \\\\  k_t=\\frac{2}{7}\\frac{m}{t_c^2}\\left(\\pi^2+(\\ln e_t)^2\\right) \\\\ c_t=-\\frac{2}{7}\\frac{m}{t_c}\\ln e_t \n\n\nwhere $k_n$, $c_n$ are normal elastic and viscous coefficients and $k_t$, $c_t$ shear elastic and viscous coefficients. For details see [Pournin2001]_.\n\n:param float m: sphere mass $m$\n:param float tc: collision time $t_c$\n:param float en: normal restitution coefficient $e_n$\n:param float es: tangential restitution coefficient $e_s$\n:return: dictionary with keys ``kn`` (the value of $k_n$), ``cn`` ($c_n$), ``kt`` ($k_t$), ``ct`` ($c_t$).");
	py::def("stressTensorOfPeriodicCell",Shop::getStress,(py::args("volume")=0),"Deprecated, use utils.getStress instead |ydeprecated|");
	py::def("normalShearStressTensors",Shop__normalShearStressTensors,(py::args("compressionPositive")=false,py::args("splitNormalTensor")=false,py::args("thresholdForce")=NaN),"Compute overall stress tensor of the periodic cell decomposed in 2 parts, one contributed by normal forces, the other by shear forces. The formulation can be found in [Thornton2000]_, eq. (3):\n\n.. math:: \\tens{\\sigma}_{ij}=\\frac{2}{V}\\sum R N \\vec{n}_i \\vec{n}_j+\\frac{2}{V}\\sum R T \\vec{n}_i\\vec{t}_j\n\nwhere $V$ is the cell volume, $R$ is \"contact radius\" (in our implementation, current distance between particle centroids), $\\vec{n}$ is the normal vector, $\\vec{t}$ is a vector perpendicular to $\\vec{n}$, $N$ and $T$ are norms of normal and shear forces.\n\n:param bool splitNormalTensor: if true the function returns normal stress tensor split into two parts according to the two subnetworks of strong an weak forces.\n\n:param Real thresholdForce: threshold value according to which the normal stress tensor can be split (e.g. a zero value would make distinction between tensile and compressive forces).");
	py::def("fabricTensor",Shop__fabricTensor,(py::args("cutoff")=0.0,py::args("splitTensor")=false,py::args("thresholdForce")=NaN),"Computes the fabric tensor $F_{ij}=\\frac{1}{n_c}\\sum_c n_i n_j$ [Satake1982]_, for all interactions $c$.\n\n:param Real cutoff: intended to disregard boundary effects: to define in [0;1] to focus on the interactions located in the centered inner (1-cutoff)^3*$V$ part of the spherical packing $V$.\n\n:param bool splitTensor: split the fabric tensor into two parts related to the strong (greatest compressive normal forces) and weak contact forces respectively.\n\n:param Real thresholdForce: if the fabric tensor is split into two parts, a threshold value can be specified otherwise the mean contact force is considered by default. Use negative signed values for compressive states. To note that this value could be set to zero if one wanted to make distinction between compressive and tensile forces.");
	py::def("bodyStressTensors",Shop__getStressLWForEachBody,"Compute and return a table with per-particle stress tensors. Each tensor represents the average stress in one particle, obtained from the contour integral of applied load as detailed below. This definition is considering each sphere as a continuum. It can be considered exact in the context of spheres at static equilibrium, interacting at contact points with negligible volume changes of the solid phase (this last assumption is not restricting possible deformations and volume changes at the packing scale).\n\nProof: \n\nFirst, we remark the identity:  $\\sigma_{ij}=\\delta_{ik}\\sigma_{kj}=x_{i,k}\\sigma_{kj}=(x_{i}\\sigma_{kj})_{,k}-x_{i}\\sigma_{kj,k}$.\n\nAt equilibrium, the divergence of stress is null: $\\sigma_{kj,k}=\\vec{0}$. Consequently, after divergence theorem: $\\frac{1}{V}\\int_V \\sigma_{ij}dV = \\frac{1}{V}\\int_V (x_{i}\\sigma_{kj})_{,k}dV = \\frac{1}{V}\\int_{\\partial V}x_i\\sigma_{kj}n_kdS = \\frac{1}{V}\\sum_bx_i^bf_j^b$.\n\nThe last equality is implicitely based on the representation of external loads as Dirac distributions whose zeros are the so-called *contact points*: 0-sized surfaces on which the *contact forces* are applied, located at $x_i$ in the deformed configuration.\n\nA weighted average of per-body stresses will give the average stress inside the solid phase. There is a simple relation between the stress inside the solid phase and the stress in an equivalent continuum in the absence of fluid pressure. For porosity $n$, the relation reads: $\\sigma_{ij}^{equ.}=(1-n)\\sigma_{ij}^{solid}$.\n\nThis last relation may not be very useful if porosity is not homogeneous. If it happens, one can define the equivalent bulk stress a the particles scale by assigning a volume to each particle. This volume can be obtained from :yref:`TesselationWrapper` (see e.g. [Catalano2014a]_)");
  py::def("getDynamicStress",Shop::getDynamicStress,  "Compute the dynamic stress tensor for each body: $\\sigma^p_D= - \\frac{1}{V^p} m^p u'^p \\otimes u'^p$");
  py::def("getTotalDynamicStress",Shop::getTotalDynamicStress, (py::args("volume")=0), "Compute the total dynamic stress tensor : $\\sigma_D= - \\frac{1}{V}\\sum_p m^p u'^p \\otimes u'^p$. The volume have to be provided for non-periodic simulations. It is computed from cell volume for periodic simulations.");
	py::def("getStress",Shop::getStress,(py::args("volume")=0),"Compute and return Love-Weber stress tensor:\n\n $\\sigma_{ij}=\\frac{1}{V}\\sum_b f_i^b l_j^b$, where the sum is over all interactions, with $f$ the contact force and $l$ the branch vector (joining centers of the bodies). Stress is negativ for repulsive contact forces, i.e. compression. $V$ can be passed to the function. If it is not, it will be equal to the volume of the cell in periodic cases, or to the one deduced from utils.aabbDim() in non-periodic cases.");
	py::def("getStressProfile",Shop::getStressProfile,(py::args("volume"),py::args("nCell"),py::args("dz"),py::args("zRef"),py::args("vPartAverageX"),py::args("vPartAverageY"),py::args("vPartAverageZ")),"Compute and return the stress tensor depth profile, including the contribution from Love-Weber stress tensor and the dynamic stress tensor taking into account the effect of particles inertia. For each defined cell z, the stress tensor reads: \n\n $\\sigma_{ij}^z= \\frac{1}{V}\\sum_c f_i^c l_j^{c,z} - \\frac{1}{V}\\sum_p m^p u'^p_i u'^p_j$,\n\n where the first sum is made over the contacts which are contained or cross the cell z, f^c is the contact force from particle 1 to particle 2, and l^{c,z} is the part of the branch vector from particle 2 to particle 1, contained in the cell. The second sum is made over the particles, and u'^p is the velocity fluctuations of the particle p with respect to the spatial averaged particle velocity at this point (given as input parameters). The expression of the stress tensor is the same as the one given in getStress plus the inertial contribution. Apart from that, the main difference with getStress stands in the fact that it gives a depth profile of stress tensor, i.e. from the reference horizontal plane at elevation zRef (input parameter) until the plane of elevation zRef+nCell*dz (input parameters), it is computing the stress tensor for each cell of height dz. For the love-Weber stress contribution, the branch vector taken into account in the calculations is only the part of the branch vector contained in the cell considered.\n To validate the formulation, it has been checked that activating only the Love-Weber stress tensor, and suming all the contributions at the different altitude, we recover the same stress tensor as when using getStress. For my own use, I have troubles with strong overlap between fixed object, so that I made a condition to exclude the contribution to the stress tensor of the fixed objects, this can be desactivated easily if needed (and should be desactivated for the comparison with getStress)." );
	py::def("getStressProfile_contact",Shop::getStressProfile_contact,(py::args("volume"),py::args("nCell"),py::args("dz"),py::args("zRef")),"same as getStressProfile, only contact contribution." );
  py::def("getDepthProfiles",Shop::getDepthProfiles,(py::args("volume"),py::args("nCell"),py::args("dz"),py::args("zRef"),py::args("activateCond"),py::args("radiusPy"),py::args("direction")),"Compute and return the particle velocity and solid volume fraction (porosity) depth profile along the direction specified (default is z; 0=>x,1=>y,2=>z). For each defined cell z, the k component of the average particle velocity reads: \n\n $<v_k>^z= \\sum_p V^p v_k^p/\\sum_p V^p$,\n\n where the sum is made over the particles contained in the cell, $v_k^p$ is the k component of the velocity associated to particle p, and $V^p$ is the part of the volume of the particle p contained inside the cell. This definition allows to smooth the averaging, and is equivalent to taking into account the center of the particles only when there is a lot of particles in each cell. As for the solid volume fraction, it is evaluated in the same way:  for each defined cell z, it reads: \n\n $<\\phi>^z= \\frac{1}{V_{cell}}\\sum_p V^p$, where $V_{cell}$ is the volume of the cell considered, and $V^p$ is the volume of particle p contained in cell z.\n This function gives depth profiles of average velocity and solid volume fraction, returning the average quantities in each cell of height dz, from the reference horizontal plane at elevation zRef (input parameter) until the plane of elevation zRef+nCell*dz (input parameters). If the argument activateCond is set to true, do the average only on particles of radius equal to radiusPy (input parameter)" );
	py::def("getDepthProfiles_center",Shop::getDepthProfiles_center,(py::args("volume"),py::args("nCell"),py::args("dz"),py::args("zRef"),py::args("activateCond"),py::args("radiusPy")),"Same as getDepthProfiles but taking into account particles as points located at the particle center." );
	py::def("getCapillaryStress",Shop::getCapillaryStress,(py::args("volume")=0,py::args("mindlin")=false),"Compute and return Love-Weber capillary stress tensor:\n\n $\\sigma^{cap}_{ij}=\\frac{1}{V}\\sum_b l_i^b f^{cap,b}_j$, where the sum is over all interactions, with $l$ the branch vector (joining centers of the bodies) and $f^{cap}$ is the capillary force. $V$ can be passed to the function. If it is not, it will be equal to one in non-periodic cases, or equal to the volume of the cell in periodic cases. Only the CapillaryPhys interaction type is supported presently. Using this function with physics MindlinCapillaryPhys needs to pass True as second argument.");
	py::def("getBodyIdsContacts",Shop__getBodyIdsContacts,(py::args("bodyID")=0),"Get a list of body-ids, which contacts the given body.");
	py::def("maxOverlapRatio",maxOverlapRatio,"Return maximum overlap ration in interactions (with :yref:`ScGeom`) of two :yref:`spheres<Sphere>`. The ratio is computed as $\\frac{u_N}{2(r_1 r_2)/r_1+r_2}$, where $u_N$ is the current overlap distance and $r_1$, $r_2$ are radii of the two spheres in contact.");
	py::def("shiftBodies",shiftBodies,(py::arg("ids"),py::arg("shift")),"Shifts bodies listed in ids without updating their velocities.");
	py::def("calm",Shop__calm,(py::arg("mask")=-1),"Set translational and rotational velocities of bodies to zero. Applied to all bodies by default. To calm only some bodies, use mask parameter, it will calm only bodies with groupMask compatible to given value");
	py::def("setNewVerticesOfFacet",setNewVerticesOfFacet,(py::arg("b"),py::arg("v1"),py::arg("v2"),py::arg("v3")),"Sets new vertices (in global coordinates) to given facet.");
	py::def("setContactFriction",Shop::setContactFriction,py::arg("angleRad"),"Modify the friction angle (in radians) inside the material classes and existing contacts. The friction for non-dynamic bodies is not modified.");
	py::def("growParticles",Shop::growParticles,(py::args("multiplier"), py::args("updateMass")=true, py::args("dynamicOnly")=true), "Change the size of spheres and clumps of spheres by the multiplier. If updateMass=True, then the mass and inertia are updated. dynamicOnly=True will select dynamic bodies.");
	py::def("growParticle",Shop::growParticle,(py::args("bodyID"),py::args("multiplier"), py::args("updateMass")=true), "Change the size of a single sphere (to be implemented: single clump). If updateMass=True, then the mass is updated.");
	py::def("intrsOfEachBody",intrsOfEachBody,"returns list of lists of interactions of each body");
	py::def("numIntrsOfEachBody",numIntrsOfEachBody,"returns list of number of interactions of each body");
	py::def("TetrahedronSignedVolume",static_cast<Real (*)(const vector<Vector3r>&)>(&TetrahedronSignedVolume),"TODO");
	py::def("TetrahedronVolume",static_cast<Real (*)(const vector<Vector3r>&)>(&TetrahedronVolume),"TODO");
	py::def("TetrahedronInertiaTensor",TetrahedronInertiaTensor,"TODO");
	py::def("TetrahedronCentralInertiaTensor",TetrahedronCentralInertiaTensor,"TODO");
	py::def("TetrahedronWithLocalAxesPrincipal",TetrahedronWithLocalAxesPrincipal,"TODO");
	py::def("momentum",Shop::momentum,"TODO");
	py::def("angularMomentum",Shop::angularMomentum,(py::args("origin")=Vector3r(Vector3r::Zero())),"TODO");
	py::def("getSpheresVolume2D",Shop__getSpheresVolume2D,(py::arg("mask")=-1),"Compute the total volume of discs in the simulation, mask parameter is considered");
	py::def("voidratio2D",Shop__getVoidRatio2D,(py::arg("zlen")=1),"Compute 2D packing void ratio $\\frac{V-V_s}{V_s}$ where $V$ is overall volume and $V_s$ is volume of disks.\n\n:param float zlen: length in the third direction.\n");
	py::def("getStressAndTangent",Shop__getStressAndTangent,(py::args("volume")=0,py::args("symmetry")=true),"Compute overall stress of periodic cell using the same equation as function getStress. In addition, the tangent operator is calculated using the equation published in [Kruyt and Rothenburg1998]_:\n\n.. math:: S_{ijkl}=\\frac{1}{V}\\sum_{c}(k_n n_i l_j n_k l_l + k_t t_i l_j t_k l_l)\n\n:param float volume: same as in function getStress\n:param bool symmetry: make the tensors symmetric.\n\n:return: macroscopic stress tensor and tangent operator as py::tuple");
	py::def("setBodyPosition",setBodyPosition,(py::args("id"),py::args("pos"),py::args("axis")="xyz"),"Set a body position from its id and a new vector3r.\n\n:param int id: the body id.\n:param Vector3 pos: the desired updated position.\n:param str axis: the axis along which the position has to be updated (ex: if axis==\"xy\" and pos==Vector3r(r0,r1,r2), r2 will be ignored and the position along z will not be updated).");
	py::def("setBodyVelocity",setBodyVelocity,(py::args("id"),py::args("vel"),py::args("axis")="xyz"),"Set a body velocity from its id and a new vector3r.\n\n:param int id: the body id.\n:param Vector3 vel: the desired updated velocity.\n:param str axis: the axis along which the velocity has to be updated (ex: if axis==\"xy\" and vel==Vector3r(r0,r1,r2), r2 will be ignored and the velocity along z will not be updated).");
	py::def("setBodyOrientation",setBodyOrientation,(py::args("id"),py::args("ori")),"Set a body orientation from its id and a new Quaternionr.\n\n:param int id: the body id.\n:param Quaternion ori: the desired updated orientation.");
	py::def("setBodyAngularVelocity",setBodyAngularVelocity,(py::args("id"),py::args("angVel")),"Set a body angular velocity from its id and a new Vector3r.\n\n:param int id: the body id.\n:param Vector3 angVel: the desired updated angular velocity.");
	py::def("setBodyColor",setBodyColor,(py::args("id"),py::args("color")),"Set a body color from its id and a new Vector3r.\n\n:param int id: the body id.\n:param Vector3 color: the desired updated color.");
}
