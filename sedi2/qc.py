#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dkunhappan <dkunhappan@dt-com012>
#
# Distributed under terms of the MIT license.

"""

"""
from __future__ import division 
from math import pi 


rho_f = 1000; 
rho_p = 7500; 
sp_rad = 2.5e-04 
d_rad = 0.007875; 
nu_f = 1e-04;
g = 10.50

sr3 = sp_rad**3 
nu2 = nu_f*nu_f 
sp_vol = (4/3)*pi*(sp_rad)**3 
d_vol = (4/3)*pi*(d_rad)**3 
num = 3000 
vol_frac = (num*sp_vol)/d_vol
eps = sp_rad/d_rad; eps3 =eps*eps*eps; 
rho_dif  = rho_p - rho_f 
Re_p = (2/9)*((rho_dif*sr3*g)/(rho_f*nu2))
Re_d = (9/2)*(vol_frac/eps3)*Re_p 

print " ...In vals.. \n  " 
print "rho_f = " , rho_f 
print "rho_p = ",  rho_p 
print "sp_rad = ", sp_rad 
print "d_rad = ", d_rad 


print "Drop Reynolds number = " , Re_d 
print "Particle Reynolds Number = ",  Re_p 
print "vol_frac =" , vol_frac 
print "epsilon = ",  eps; 
