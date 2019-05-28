# Yade-OpenFOAM-coupling
An OpenFOAM solver for realizing CFD-DEM simulations with the Open Source Discrete Element Solver Yade-DEM. 
 * Fast mesh search (including range based search) based on k-d Tree algorithm, faster than the original Octree search offered by OpenFOAM (mesh.findCell,  mesh.findNearest).
 * Gaussian interpolation of field variables. 
 * Simple point-force coupling (icoFoamYade) solver
 * Full 4-way coupling (pimpleFoamYade) solver.
 * Documentation : https://yade-dem.org/doc/FoamCoupling.html
 * Examples : https://gitlab.com/yade-dev/trunk/tree/master/examples/openfoam


![Alt text](ccpl1.png)

Prerequisites : Latest Yade git version with the FoamCoupling engine (https://gitlab.com/yade-dev/trunk). OpenFOAM-6. 

## Build
* Compile the lib FoamYade and solver icoFoamYade : 
  * ``./Allwmake``

## Running 
* Copy example to $FOAM_RUN
* Create a symbolic link
  * ``ln -s /path/to/your/yade/install/bin/yade-exec libyade.py``
* Run 
  * ``cp -r example_icoFoamYade /to/your/run/dir ``
  * ``mpiexec -n 1 python scriptYade.py : -n 2 icoFoamYade -parallel``
  
