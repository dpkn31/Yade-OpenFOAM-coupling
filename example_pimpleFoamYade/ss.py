numspheres = 1000; 
young = 5e6; 
density =1000; 

O.materials.append(FrictMat(young=young,poisson=0.5,frictionAngle=radians(15),density=density,label='spheremat'))
O.materials.append(FrictMat(young=young*100,poisson=0.5,frictionAngle=radians(15),density=0,label='wallmat'))

val = 2e-08

f1 = facet([Vector3(0,val,0), Vector3(1,val,0), Vector3(0,val,1)], material='wallmat'); 
f2 = facet([Vector3(0,val,1), Vector3(1,val,1), Vector3(1,val,0)], material='wallmat'); 
O.bodies.append(f1); O.bodies.append(f2); 

#f3 = facet([Vector3(0,1-val,0), Vector3(1,1-val,0), Vector3(0,1-val,1)], material='wallmat'); 
#f4 = facet([Vector3(0,1-val,1), Vector3(1,1-val,1), Vector3(1,1-val,0)], material='wallmat'); 
#O.bodies.append(f3); O.bodies.append(f4); 

sd = Vector3(0.5,1.5,0.5); rad=0.15; 

#sp = pack.SpherePack();
sp = pack.randomDensePack(pack.inSphere(sd, rad), spheresInCell=60, radius=0.01, rRelFuzz = 0.0, returnSpherePack=True)

O.bodies.append([sphere(center,rad,material='spheremat') for center,rad in sp]) 
idlist = [b.id for b in O.bodies if b.shape==Sphere]

newton=NewtonIntegrator(damping=0.0, gravity = (0.0 ,-9.81, 0.0)) # add small damping in case of stability issues.. ~ 0.1 max, also note : If using gravity,  make sure buoyancy force is switched on in the foam side, else the simulation will crash


O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb()], allowBiggerThanPeriod=True),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(),Ig2_Facet_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_FrictPhys()],
		[Law2_ScGeom_FrictPhys_CundallStrack()]
	),
	GlobalStiffnessTimeStepper(timestepSafetyCoefficient=0.1, label = "ts"), 
	newton, 
	#VTKRecorder(fileName='yadep/3d-vtk-',recorders=['spheres'],iterPeriod=1000)
]

