#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

from __future__ import absolute_import, division, print_function
import numpy as np

def stv_handl_spring(imode, statev, eelas, eplas, eqplas):
    """
    This function stores (imode = 0) and returns (imode = 1) elastic deformation,
    plastic deformation, and Alpha variable (internal hardening variable)
    
    """    
    if imode == 0:
        eelas  = statev[0]
        eplas  = statev[1]
        eqplas = statev[2]
    else:
        statev[0] = eelas
        statev[1] = eplas
        statev[2] = eqplas
    
    return statev, eelas, eplas, eqplas

def umat_spring(stress, strann, statev, par):
    """   
    This function contents the formulation of a nonlinear formulation of
    a 1D frictional sring with isotropic hardening. 
    
    Algorithm of 1-D Rate independent plasticity, Isotropic hardening.
    Chapter 1. Box 1.4, Simo & Hughes. Computational Inelasticity.
    
    Parameters
    ----------
    statev: ndarray
            State variables array at the current integration point. 
    stress: ndarray
            Stress tensor at the current integration point.
    strann: nd array
            Strain tensor at the current integration point.
    par   : nd array
            Material properties or parameters of the element i

    Returns
    -------
    E     : Modulus of elasticity updated, considering the isotropic hardening plastic model.
    
    Other functions associated:
    --------------------------
    stv_handl_spring : Stores and returns state variables of the material elastoplastic plastic model behavior

    """    
    
    Emod       = par[2]
    kh         = par[7]
    sigy       = par[8]
    
    eelas   = 0.0
    eplas   = 0.0
    eqplas  = 0.0

    statev, eelas, eplas, eqplas = stv_handl_spring(0, statev, eelas, eplas, eqplas)
        
    sig_trial = Emod*(strann - eplas)
    ftrial = np.abs(sig_trial)-(sigy+kh*eqplas)
        
    if ftrial > 0.0:
        gama_par = (ftrial/(Emod+kh))
        stress   = (1.0 -(gama_par*Emod/sig_trial))*sig_trial
        eplas    = eplas + gama_par*np.sign(sig_trial)
        eqplas   = eqplas + gama_par
        eelas    = strann - eplas
        E = Emod*kh/(Emod+kh)
    else:
        stress = sig_trial
        eplas  = eplas
        eqplas = eqplas
        eelas  = eelas
        E = Emod

    statev, eelas, eplas, eqplas = stv_handl_spring(1, statev, eelas, eplas, eqplas)   
    return stress , statev , E


def stv_handl_PileRCK(imode, statev, Ye, Yp):
    """
    This function stores (imode = 0) and returns (imode = 1) elastic deformation Ye and
    plastic deformation Yp
    
    Contistitive model of a P-y curve for stiff soil (rock).
    Handbook of design of piles under later load. Section 2.9
    
    """    
    if imode == 0:
        Ye  = statev[0]
        Yp  = statev[1]
    else:
        statev[0] = Ye
        statev[1] = Yp

    return statev, Ye, Yp

def umat_PileRCK(Pi, Yi, statev, par):
    """   
    This function contents the formulation of a nonlinear formulation of
    a 1D spring in X direction with P-y curves.
    
    Contistitive model of a P-y curve for stiff soil (rock).
    Handbook of design of piles under later load. Section 2.9
    
    Parameters
    ----------
    Pi:     float
            Load at current time/deformation increment.
    Yi:     float
            Displacement at current time/deformation increment.
    statev: ndarray
            State variables array at the current increment.        
    par:    ndarray
            Material properties or parameters of the element i

    Returns
    -------
    E     : Modulus of elasticity updated, considering the isotropic hardening plastic model.
    
    Other functions associated:
    --------------------------
    stv_handl_spring : Stores and returns state variables of the material elastoplastic plastic model behavior

    """    
    
    K1   = par[0]
    Yp0  = par[1]
    K2   = par[2]*K1
    Pult = par[3]
    Py   = K1*Yp0
    
    Ye = 0.0
    Yp = 0.0
    statev, Ye, Yp = stv_handl_PileRCK(0, statev, Ye, Yp)
    
    if Yi > Yp0:
       E = K2
       P = Py + E*(Yi-Yp0)
       #
       if (Py + P) >= Pult:
          if K2 == 0.0:
             P = Py
          else:
             print ('Spring has reached maximum capacity')  
             P = Pult
          #End if   
          E = 0.0
       #End if   
    else:
       Ye = Yi
       E = K1
       P = E*Yi
    # End if
    statev, Ye, Yp = stv_handl_PileRCK(1, statev, Ye, Yp)    

    return P, statev, E


def stv_handl_PileCLY(imode, statev, Pant, Yant):
    """
    This function stores (imode = 0) and returns (imode = 1) elastic deformation Ye and
    plastic deformation Yp
    
    Contistitive model of a P-y curve for clay.
    Handbook of design of piles under later load. Section 2.4
    Matlock, 1970. P-y curve
    
    """    
    if imode == 0:
        Pant  = statev[0]
        Yant  = statev[1]
    else:
        statev[0] = Pant
        statev[1] = Yant

    return statev, Pant, Yant

def umat_PileCLY(Pi, Yi, statev, par):
    """   
    This function contents the formulation of a nonlinear formulation of
    a 1D spring in X direction with P-y curves.
    
    Contistitive model of a P-y curve for clay soil.
    Handbook of design of piles under later load. Section 2.4
    Matlock, 1970. P-y curve
    
    Parameters
    ----------
    Pi:     float
            Load at current time/deformation increment.
    Yi:     float
            Displacement at current time/deformation increment.
    statev: ndarray
            State variables array at the current increment.        
    par:    ndarray
            Material properties or parameters of the element i

    Returns
    -------
    E     : Modulus of elasticity updated, considering the isotropic hardening plastic model.
    
    Other functions associated:
    --------------------------
    stv_handl_spring : Stores and returns state variables of the material elastoplastic plastic model behavior

    """    
    #
    c     = par[0]
    b     = par[1]
    Gamma = par[2]
    x     = par[3]
    J     = par[4]
    e50   = par[5]
    #
    Pu1   = 9*c*b
    Pu2   = (3 + (Gamma*x/c) + (J*x/b))*c*b
    Pu    = min(Pu1, Pu2)
    Y50   = 2.5*e50*b
    Ymax  = 8*Y50
    #  
    Pant  = 0.0
    Yant  = 0.0
    statev, Pant, Yant = stv_handl_PileCLY(0, statev, Pant, Yant)
    #
    if Yi > Y50:
       P = 0.5*Pu*(Yi-Y50)**(1/3)
       #
       if Yi >= Ymax:
          print ('Spring has reached maximum capacity')
          P = Pant
       #End if   
    else:
       Eelas = Pu*0.5/Y50
       P = Eelas*Yi
    # End if
    #
    E    = (P-Pant) / (Yi-Yant)
    Pant = P
    Yant = Yi
    #
    statev, Pant, Yant = stv_handl_PileCLY(1, statev, Pant, Yant)    
    #
    return P, statev, E
























