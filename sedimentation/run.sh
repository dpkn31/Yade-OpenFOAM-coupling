#! /bin/bash
mpiexec -n 1 python scriptYade.py : -n 4 pimpleFoamYade -parallel

