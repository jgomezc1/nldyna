# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 15:28:36 2019

@author: AX201 GMRS
"""

from __future__ import absolute_import, division, print_function
import numpy as np
import uel_solid as uel
#
# Problem Parameters
#
nprops = 7
nsvars = 88
ntens  = 4
props   = np.zeros([nprops])
svars   = np.zeros([nsvars])
du = np.zeros([8])

coord =([0.0 , 0.0], [1.0 , 0.0], [1.0 , 1.0], [0.0 , 1.0])

props[0]= 52.0e3
props[1]= 0.33
props[2]= 60.0
props[3]= 37.0
props[4]= 383.3
props[5]= 2040.0
props[6]= 1000.0

kl, ml, cl, svars, rhsl = uel.uel4nquad_PLK(coord, props , svars , du)
print(kl)
