#! /bin/bash
mpiexec -n 1 python scriptYade.py : -n 2 pimpleFoamYade -parallel

