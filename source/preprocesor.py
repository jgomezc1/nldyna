# -*- coding: utf-8 -*-
"""
This module contains functions to preprocess the input files to compute
a Finite Element Analysis.

"""
from __future__ import division, print_function
import numpy as np
import femutil as fem
        
def readin(folder):
    """Read the input files"""
    inipar        = np.loadtxt(folder + '01_Inipar.txt', ndmin=2, usecols=(0), skiprows=6)
    nodes         = np.loadtxt(folder + '03_Nodes.txt' , ndmin=2, skiprows=3)
    loads         = np.loadtxt(folder + '05_Nodal_loads.txt' , ndmin=2, skiprows=3)
    #
    # -------------------------------------------------------------------------------------
    # Read constraints file and identify if there is any diaphfragm or constraint
    const         = np.loadtxt(folder + '07_DOF_Constraints.txt' , ndmin=2, skiprows=9)
    DPH_flag = 0
    CST_flag = 0
    #
    if len(const) != 0:
       if max(const[:,1]) >= 0:
           DPH_flag = 1
       if min(const[:,1]) < 0:
           CST_flag = 1
       # End if
    #End if   
    const = list(const)
    const.insert(0,CST_flag)
    const.insert(0,DPH_flag)
    #
    # -------------------------------------------------------------------------------------
    #
    NLSTA  = int(inipar[5,0])
    NLDYNA = int(inipar[6,0])
    #
    if (NLSTA == 1) and (NLDYNA == 0):
       Seismo_signal = []
    elif (NLSTA == 0) and (NLDYNA == 1):
       Seismo_signal = np.loadtxt(folder + '06_Seismo_signal.txt' , skiprows=5)
    elif (NLSTA == 1) and (NLDYNA == 1):
       Seismo_signal = np.loadtxt(folder + '06_Seismo_signal.txt' , skiprows=5)
    # End if
    #
    # -------------------------------------------------------------------------------------   
    # Read material - parameters file
    mat_file = open(folder + '02_Mater.txt', 'r')
    #
    mats    = []
    matflag = 0
    for line in mat_file:
         if matflag == 1:
              l = np.array(line.split(), dtype = float)
              mats.append(l)
         # End if
         if line == 'PARAMETERS INFORMATION\n':
            matflag = 1
         #End if
    #End for
    #
    # -------------------------------------------------------------------------------------        
    # Read elements file
    eles_file = open(folder + '04_Eles.txt', 'r')
    #
    elements = []
    eleflag  = 0
    for line in eles_file:
         if eleflag == 1:
              l = np.array(line.split(), dtype = int)
              elements.append(l)
         #End if
         if line == 'ELEMENTS\n':
            eleflag = 1
         #End if
    #End for
    #
    # -------------------------------------------------------------------------------------        
    # Initialize Msvar (List where state variables of each element will be stocked)
    Msvar = []
    ILF   = []
    nele  = len(elements)
    #
    for i in range (nele):
       #
       iet = elements[i][1]
       ndof, nnodes, ngpts = fem.eletype(iet)
       #
       if iet == 0: # Linear - 1D Spring
          nsvar = 1
       if iet == 1: # Linear - Simple frame
          nsvar = 1
       if iet == 2: # Linear - Full frame
          nsvar = 1
       if iet == 3: # Linear - 2D Truss
          nsvar = 1
       if iet == 4: # NonLin - 4 noded plate (Rutinas de Jefe)
          nsvar = 88
       if iet == 5: # NonLin - 1D Spring
          nsvar = 5
       if iet == 6: # Lin - 2D shear-rotational Spring
          nsvar = 1
       if iet == 7: # Lin - 1D rotational Spring
          nsvar = 1 
       if iet == 8: # NonLin - 1D rotational Spring
          nsvar = 5 
       if iet == 9: # NonLin - 1D pile spring for stiff soils
          nsvar = 4
       if iet == 10: # Linear - 3D full frame bending ans shear effects
          nsvar = 1    
       # End if
       #
       ele_svar = np.zeros((nsvar))
       ele_ilf  = np.zeros((ndof))
       Msvar.append(ele_svar)
       ILF.append(ele_ilf)
    # End for
    # -------------------------------------------------------------------------------------   
    return inipar, nodes, mats, elements, loads, Msvar, ILF, Seismo_signal, const

def intparams(inipar):
    #
    ac = np.zeros(8, dtype = float)
    dt = inipar[0,0]
    T  = inipar[1,0]
    #
    m     = int(T/dt)
    theta = 1.42
    #
    ac[0] = 6.0/(theta*dt)
    ac[1] = theta*dt/2.0
    ac[2] = 3.0/(theta*dt)
    ac[3] = 6.0/(theta*theta*dt*dt)
    ac[4] = 1.0/theta
    ac[5] = dt/2.0
    ac[6] = dt*dt/2.0
    ac[7] = dt*dt/6.0
    #
    return  m, T, dt, ac, theta

