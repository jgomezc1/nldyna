# -*- coding: utf-8 -*-
"""
Created on Wed May 15 06:48:28 2019

@author: JULIAN PARRA
"""

import numpy as np

def eigen(inipar, M, K):
    """
    This function computes eigenvalues and eigenvectors of the system assembly
    
    Parameters:
    -----------
    inipar : Global parameters
    M      : Mass matrix     
    K      : Stiffness matrix
    
    """
    NLSTA  = int(inipar[5,0])
    NLDYNA = int(inipar[6,0])
    
    if (NLSTA == 1) and (NLDYNA == 0):
         T = 'Not computed,static system solution'
         Tmin = ''
         #
    elif (NLSTA == 0) and (NLDYNA == 1):
         Minv = np.linalg.inv(M)
         w2   = np.dot(Minv,K)    ### Sin llevar a valores ppales, frecuencias al cuadrado
         eigenvalues, eigenvectors  = np.linalg.eig(w2) ### Valores y vectores propios
         w = eigenvalues**0.5
         T = 2*np.pi/w
         T[::-1].sort()
         Tmin = min(T)
         #
    elif (NLSTA == 1) and (NLDYNA == 1):
         Minv = np.linalg.inv(M)
         w2   = np.dot(Minv,K)    ### Sin llevar a valores ppales, frecuencias al cuadrado
         eigenvalues, eigenvectors  = np.linalg.eig(w2) ### Valores y vectores propios
         w = eigenvalues**0.5
         T = 2*np.pi/w
         T[::-1].sort()
         Tmin = min(T)
         #
    # End if  
    #
    return T, Tmin
     