import sys 
from libyade import *
from yade.utils import *
from numpy import mod
initMPI() 


class simulation(): 
   
  def __init__(self):
    
    
    
    numspheres = 20000; 
    young = 5e6; density = 1000;

    
    O.materials.append(FrictMat(young=1e6,poisson=0.5,frictionAngle=radians(15),density=1000,label='spheres'))
    O.materials.append(FrictMat(young=1e6,poisson=0.5,frictionAngle=0,density=0,label='walls'))
    
    mn, mx= Vector3(0,0,0), Vector3(0.999,0.999, 0.999)
    walls=aabbWalls([mn,mx],thickness=0,material='walls')
    wallIds=O.bodies.append(walls)
 
    sp = pack.SpherePack();
    sp.makeCloud(mn,mx,-1,0.3333,numspheres,False, 0.95,seed=1) 
    O.bodies.append([sphere(center,rad,material='spheres') for center,rad in sp]) 

    sphereIDs = [b.id for b in O.bodies if type(b.shape)==Sphere] 
    utils.setnumParticles(numspheres)
    utils.setBodylist(sphereIDs) 

    self.dt_yade = 8.00e-05; self.dt_foam = 0.05; 
    self.dataex = int (self.dt_foam/self.dt_yade); 
    
     

    newton=NewtonIntegrator(damping=0)

    O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_FrictPhys()],
		[Law2_ScGeom_FrictPhys_CundallStrack()]
	),
        PyRunner(command='sim.coupling()', iterPeriod=1, label='coupling'),
        PyRunner(command='sim.sethydroforce()', iterPeriod=1, label='hydroforce'), 
	#GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.7), #TODO:time coupling based on stokes num. 
	newton, 
        VTKRecorder(fileName='yadep/3d-vtk-',recorders=['spheres'],iterPeriod=self.dataex*10)
    ]



  def sethydroforce(self):
     if (O.iter==1 or (mod(O.iter, 1000)==0)):
       self.printMessage();  
     setHydroForce(); 

    
        
  def coupling(self): 
      if (O.iter==1 or (mod(O.iter, self.dataex)==0)):
       runCoupling()
      
         

  def irun(self,num):
      O.dt = self.dt_yade 
      O.run(num,1)


  def killmpi(self):
      killMPI(); 



  def printMessage(self):
     print "********************************YADE-ITER = " + str(O.iter) +" **********************************" 
   




if __name__=="__main__":
    sim = simulation()
    sim.irun(10000000)
    sim.killMPI()

import __builtin__ 
__builtin__.sim=sim
