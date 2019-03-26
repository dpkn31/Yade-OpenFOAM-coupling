#Deepak Kunhappan, deepak.kunhappan@3sr-grenoble.fr 
#Example script of Yade-OpenFOAM couplinfg.
#get the OpenFOAM solver at : https://github.com/dpkn31/Yade-OpenFOAM-coupling
#get the latest version of Yade with the FoamCoupling engine here : https://gitlab.com/yade-dev/trunk


import sys 
from libyade import *
from yade.utils import *

initMPI() #Initialize the mpi environment, always required. 
fluidCoupling = FoamCoupling();
fluidCoupling.getRank();


#example of spheres in shear flow : two-way point force coupling
class simulation(): 
   
  def __init__(self):

    O.periodic = True; 
    O.cell.setBox(0.012,0.012,0.012); 
    
    
    young = 5e6; density = 3000;

    O.materials.append(FrictMat(young=young,poisson=0.5,frictionAngle=radians(15),density=density,label='spheremat'))
    O.materials.append(FrictMat(young=young,poisson=0.5,frictionAngle=0,density=0,label='wallmat'))
   
    
    
    sd = Vector3(0.006,0.01,0.006); rad=1.2e-03;
 
    sp = pack.randomDensePack(pack.inSphere(sd, rad), spheresInCell=50, radius=0.60e-04, rRelFuzz = 0.0,returnSpherePack=True)
    O.bodies.append([sphere(center,rad,material='spheremat') for center,rad in sp]) 

    sphereIDs = [b.id for b in O.bodies if type(b.shape)==Sphere] 
    numparts = len(sphereIDs);     


    fluidCoupling.setNumParticles(numparts)
    fluidCoupling.setIdList(sphereIDs)
    fluidCoupling.isGaussianInterp=True;  #use pimpleFoamYade for gaussianInterp

    newton=NewtonIntegrator(damping=0.0, gravity = (0.0 ,0.0, 0.0)) # add small damping in case of stability issues.. ~ 0.1 max, also note : If using gravity,  make sure buoyancy force is switched on in the foam side, else the simulation will crash

    O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_FrictPhys()],
		[Law2_ScGeom_FrictPhys_CundallStrack()]
	),
	GlobalStiffnessTimeStepper(timestepSafetyCoefficient=0.7, label = "ts"),
        VTKRecorder(fileName='yadep/3d-vtk-',recorders=['spheres'],iterPeriod=1000), 
        fluidCoupling, #to be called after timestepper  
        PyRunner(command='sim.printMessage()', iterPeriod= 10, label='hydroforce'), 
	newton, 

    ]

  def printMessage(self):
     print "********************************YADE-ITER = " + str(O.iter) +" **********************************" 



  def irun(self,num): 
      O.run(num,1)


if __name__=="__main__":
    sim = simulation()
    sim.irun(1000000)
    fluidCoupling.killmpi()

import __builtin__ 
__builtin__.sim=sim
