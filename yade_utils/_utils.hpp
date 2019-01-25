#pragma once

#include<pkg/dem/Shop.hpp>
#include<core/Scene.hpp>
#include<core/Omega.hpp>
#include<pkg/dem/ScGeom.hpp>
#include<pkg/dem/DemXDofGeom.hpp>
#include<pkg/common/Facet.hpp>
#include<pkg/dem/Tetra.hpp>
#include<pkg/common/Sphere.hpp>
#include<pkg/common/NormShearPhys.hpp>
#include<lib/computational-geometry/Hull2d.hpp>
#include<lib/pyutil/doc_opts.hpp>
#include<pkg/dem/ViscoelasticPM.hpp>
#include <mpi.h> 
#include <vector>

namespace py = boost::python;


void initMPI(); 
void setnumParticles(int);
void setBodylist(const py::list& );
void castParticle(); 
void updateProclist(); 
void resetProcList(); 
void recvHydroForce();
void killMPI(); 
void runCoupling();
void setHydroForce();
void castNumpParticle(int& ); 

py::tuple negPosExtremeIds(int axis, Real distFactor=1.1);

py::tuple coordsAndDisplacements(int axis,py::tuple Aabb=py::tuple());

void setRefSe3();

Real PWaveTimeStep();
Real RayleighWaveTimeStep();

py::tuple interactionAnglesHistogram(int axis, int mask=0, size_t bins=20, py::tuple aabb=py::tuple(), bool sphSph=0, Real minProjLen=1e-6);

py::tuple bodyNumInteractionsHistogram(py::tuple aabb=py::tuple());

Vector3r inscribedCircleCenter(const Vector3r& v0, const Vector3r& v1, const Vector3r& v2);
py::dict getViscoelasticFromSpheresInteraction(Real tc, Real en, Real es);

/* reset highlight of all bodies */
void highlightNone();

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
Real sumTorques(py::list ids, const Vector3r& axis, const Vector3r& axisPt);

/* Sum forces acting on bodies within mask.
 *
 * @param ids list of ids
 * @param direction direction in which forces are summed
 *
 */
Real sumForces(py::list ids, const Vector3r& direction);

/* Sum force acting on facets given by their ids in the sense of their respective normals.
   If axis is given, it will sum forces perpendicular to given axis only (not the the facet normals).
*/
Real sumFacetNormalForces(vector<Body::id_t> ids, int axis=-1);

/* Set wire display of all/some/none bodies depending on the filter. */
void wireSome(string filter);
void wireAll();
void wireNone();
void wireNoSpheres();


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
bool pointInsidePolygon(py::tuple xy, py::object vertices);

/* Compute area of convex hull when when taking (swept) spheres crossing the plane at coord, perpendicular to axis.

	All spheres that touch the plane are projected as hexagons on their circumference to the plane.
	Convex hull from this cloud is computed.
	The area of the hull is returned.

*/
Real approxSectionArea(Real coord, int axis);

/* Find all interactions deriving from NormShearPhys that cross plane given by a point and normal
	(the normal may not be normalized in this case, though) and sum forces (both normal and shear) on them.
	
	Returns a 3-tuple with the components along global x,y,z axes, which can be viewed as "action from lower part, towards
	upper part" (lower and upper parts with respect to the plane's normal).

	(This could be easily extended to return sum of only normal forces or only of shear forces.)
*/
Vector3r forcesOnPlane(const Vector3r& planePt, const Vector3r&  normal);

/* Less general than forcesOnPlane, computes force on plane perpendicular to axis, passing through coordinate coord. */
Vector3r forcesOnCoordPlane(Real coord, int axis);

py::tuple spiralProject(const Vector3r& pt, Real dH_dTheta, int axis=2, Real periodStart=std::numeric_limits<Real>::quiet_NaN(), Real theta0=0);

shared_ptr<Interaction> Shop__createExplicitInteraction(Body::id_t id1, Body::id_t id2);

Real Shop__unbalancedForce(bool useMaxForce /*false by default*/);
py::tuple Shop__totalForceInVolume();
Real Shop__getSpheresVolume(int mask=-1);
Real Shop__getSpheresMass(int mask=-1);
py::object Shop__kineticEnergy(bool findMaxId=false);

Real maxOverlapRatio();

Real Shop__getPorosity(Real volume=-1);
Real Shop__getVoxelPorosity(int resolution=200, Vector3r start=Vector3r(0,0,0),Vector3r end=Vector3r(0,0,0));

//Matrix3r Shop__stressTensorOfPeriodicCell(bool smallStrains=false){return Shop::stressTensorOfPeriodicCell(smallStrains);}
py::tuple Shop__fabricTensor(Real cutoff=0.0,bool splitTensor=false, Real thresholdForce=NaN);
py::tuple Shop__normalShearStressTensors(bool compressionPositive=false, bool splitNormalTensor=false, Real thresholdForce=NaN);

py::list Shop__getStressLWForEachBody();

py::list Shop__getBodyIdsContacts(Body::id_t bodyID=-1);

Real shiftBodies(py::list ids, const Vector3r& shift);

void Shop__calm(int mask=-1);

void setNewVerticesOfFacet(const shared_ptr<Body>& b, const Vector3r& v1, const Vector3r& v2, const Vector3r& v3);

py::list intrsOfEachBody();

py::list numIntrsOfEachBody();

/* The 5 following setters are used to workaround a long-standing bug in the c++/python binding which produces a memory leak (see two links below).
 * https://bugs.launchpad.net/yade/+bug/1041084
 * https://answers.launchpad.net/yade/+question/253112
 * It is not in the spirit of Yade Python binding but you can use them if you massively update bodies attributes.
 * TODO : remove them as soon as the bug is solved.
*/

/* Set a body position from its id and a new vector3r.
 *  @param id is the body id
 *  @param newPos is the desired updated position
 *  @param axis is the axis along which the position has to be updated (ex: if axis=="xy" and newPos==Vector3r(r0,r1,r2), r2 will be ignored and the position along z will not be updated).
*/
void setBodyPosition(int id, Vector3r newPos, string axis="xyz");

/* Set a body velocity from its id and a new vector3r.
 *  @param id is the body id
 *  @param newPos is the desired updated velocity
 *  @param axis is the axis along which the velocity has to be updated (ex: if axis=="xy" and newVel==Vector3r(r0,r1,r2), r2 will be ignored and the velocity along z will not be updated).
*/
void setBodyVelocity(int id, Vector3r newVel, string axis="xyz");

/* Set a body orientation from its id and a new Quaternionr.
 *  @param id is the body id
 *  @param newOri is the desired updated orientation
*/
void setBodyOrientation(int id, Quaternionr newOri);

/* Set a body angular velocity from its id and a new Vector3r.
 *  @param id is the body id
 *  @param newAngVel is the desired updated angular velocity
*/
void setBodyAngularVelocity(int id, Vector3r newAngVel);

/* Set a body color from its id and a new Vector3r.
 *  @param id is the body id
 *  @param newColor is the desired rgb color
*/
void setBodyColor(int id, Vector3r newColor);
