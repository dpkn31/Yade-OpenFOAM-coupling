# Yade-OpenFOAM-coupling
![Alt text](ccpl1.png)
A coupling module for coupling  the DEM solver YADE with the FVM solver OpenFOAM 

Prerequisites : Yade version 2019-01-04.git-f5aa5f7. OpenFOAM-6. 

* Compile icoFoamYade : 
  * icoFoamYade/FoamYade/commYade
    * Make the lib : wmake 
  * icoFoamYade/FoamYade 
    * Make the lib : wmake 
  * icoFoamYade:
    * wmake

* Copy _utils.cpp and _utils.hpp to your YADE trunk/py/ directory and then compile. 

## Running 
* Copy example to $FOAM_RUN
* Create a symbolic link
  * ln -s /path/to/your/yade/install/bin/yade-exec libyade.py
* Run 
  * mpiexec -n 1 python scriptYade.py : -n 4 icoFoamYade -parallel
  
