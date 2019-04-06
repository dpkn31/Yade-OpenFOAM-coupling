numspheres = 1000;
young = 5e6;
density =1000;

O.materials.append(FrictMat(young=young,poisson=0.5,frictionAngle=radians(15),density=density,label='spheremat'))
O.materials.append(FrictMat(young=young*100,poisson=0.5,frictionAngle=radians(15),density=0,label='wallmat'))

val = 2e-08

# f1 = facet([Vector3(0,val,0), Vector3(1,val,0), Vector3(0,val,1)], material='wallmat');
# f2 = facet([Vector3(0,val,1), Vector3(1,val,1), Vector3(1,val,0)], material='wallmat');
# O.bodies.append(f1); O.bodies.append(f2);

#f3 = facet([Vector3(0,1-val,0), Vector3(1,1-val,0), Vector3(0,1-val,1)], material='wallmat');
#f4 = facet([Vector3(0,1-val,1), Vector3(1,1-val,1), Vector3(1,1-val,0)], material='wallmat');
#O.bodies.append(f3); O.bodies.append(f4);
mn, mx = Vector3(0.04,0.16,0.04), Vector3(0.08,0.20,0.08)
sd = Vector3(0.06,0.18,0.06); rad=0.0075;
radius = 0.00025;

sp = pack.SpherePack();
sp.makeCloud(mn,mx,rMean=radius,rRelFuzz=0.00, num=48000)
#sp=pack.randomDensePack(pack.inAlignedBox(mn,mx),spheresInCell=100,radius=radius,rRelFuzz=0,returnSpherePack=True)
# sp = pack.randomDensePack(pack.inSphere(sd, rad), radius=radius, rRelFuzz = 0.0, returnSpherePack=True)
cc = [x[0] for x in sp];
coords = [];

for pos in cc:
	dx =  (pos-sd).norm()
	if dx <= 0.55*rad or dx > 1.35*rad:
		pass
	else:
		coords.append(pos)


print "len of coords = ", len(coords)
O.bodies.append([sphere(center,radius,material='spheremat') for center in coords])
idlist = [b.id for b in O.bodies if b.shape==Sphere]

fl = open('coords_r=0pt005.txt', 'w')
for j in coords:
	fl.write('%12.5e %12.5e %12.5e \n' % (j[0],j[1],j[2]))
fl.close()


newton=NewtonIntegrator(damping=0.0, gravity = (0.0 ,0.0, 0.0)) # add small damping in case of stability issues.. ~ 0.1 max, also note : If using gravity,  make sure buoyancy force is switched on in the foam side, else the simulation will crash


O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(),Ig2_Facet_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_FrictPhys()],
		[Law2_ScGeom_FrictPhys_CundallStrack()]
	),
	GlobalStiffnessTimeStepper(timestepSafetyCoefficient=0.1, label = "ts"),
	newton,
	#VTKRecorder(fileName='yadep/3d-vtk-',recorders=['spheres'],iterPeriod=1000)
]
