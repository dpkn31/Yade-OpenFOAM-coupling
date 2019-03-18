# Yade-OpenFOAM-coupling
![Alt text](ccpl1.png)
A coupling module for coupling  the DEM solver YADE with the FVM solver OpenFOAM. 

Prerequisites : Yade version 2019-01-04.git-f5aa5f7. OpenFOAM-6. 

## Build
* Compile everything : 
  * ``./Allwmake``

## Running 
* Copy example to $FOAM_RUN
* Create a symbolic link
  * ``ln -s /path/to/your/yade/install/bin/yade-exec libyade.py``
* Run 
  * ``mpiexec -n 1 python scriptYade.py : -n 4 icoFoamYade -parallel``
  
