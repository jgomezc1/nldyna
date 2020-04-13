# -*- coding: utf-8 -*-
"""
Element subroutines
-------------------
Each UEL subroutine computes the local stiffness, mass and damping matrixes for a given
finite element.

New elements can be added by including additional subroutines.

"""
from __future__ import division, print_function
import numpy as np
import umatutil as umat

#####################        NON-LINEAR 1D FRICTIONAL SPRING        ######################

def svarshandl_spring(imode , svars , statev , stress , strann):
    """Handles state variables array svars according to imode
    
    Parameters
    ----------
    imode  : int
           Handling mode:
           (0) Retrieves state variables
           (1) Updates state variables
    svars  : ndarray
           State variables array at all the integration points to be
           retrieved or updated according to imode
    statev: ndarray
          State variables array at the current integration point. 
    stress: ndarray
          Stress tensor at the current integration point.
    strann: nd array
          Strain tensor at the current integration point.
    igp   : int
          Pointer to the svars block for the current integration point.
    Returns
    -------
    svars , statev, stress, strann
    """
    
    if imode == 0:
        stress    = svars[0]
        strann    = svars[1]
        statev[0] = svars[2]
        statev[1] = svars[3]
        statev[2] = svars[4]
    else:
        svars[0]  = stress
        svars[1]  = strann
        svars[2]  = statev[0]
        svars[3]  = statev[1]
        svars[4]  = statev[2]
     # End if
    return svars, statev, stress, strann

def NL_1Dspring(iele_disp, coord, par, svar):
    """ Non-Linear 1D Friccional spring element 
    
    Algorithm of 1-D Rate independent plasticity, Isotropic hardening.
    Chapter 1. Box 1.4, Simo & Hughes. Computational Inelasticity.
    
    Parameters
    ----------
    iele_disp : ndarray
                Nodal global displacement values at time t 
    coord     : ndarray
                Coordinates for the nodes of the element (2, 2).
    par       : ndarray
                Material properties
    svar      : 1D array 
                Element's plastic variables of the increment    
    Returns
    -------
    kG   : ndarray
           Stiffness matrix in global coordinates for the element (2, 2).
    mG   : ndarray
           Mass matrix in global coordinates for the element      (2, 2).
    cG   : ndarray
           Damping matrix in global coordinates for the element   (2, 2).
    svar : 1D array 
           Element's plastic variables of the increment.
    """
    
    A      = par[0]
    Emod   = par[2]
    rho    = par[3]
    calpha = par[4]
    cbeta  = par[5]
    g      = par[6]
    vec    = coord[1, :] - coord[0, :]
    L      = np.linalg.norm(vec)
    Q      = Q_1Dspring(vec)
    #
    strann = 0.0
    stress = 0.0
    statev = np.zeros([5])
    #
    svar, statev, stress, strann = svarshandl_spring(0, svar, statev, stress, strann)
    strann = strann + (iele_disp[1] - iele_disp[0])*(1.0/L)
    stress, statev, Emod = umat.umat_spring(stress, strann, statev, par)
    #
    kl   = (A*Emod/L)*np.array([[ 1.0,-1.0],
                                [-1.0, 1.0]])
    #
    ilf  = np.dot(kl, np.dot(Q,iele_disp))
    kG   = np.dot(np.dot(Q.T, kl), Q)
    #
    mG   = (0.5*A*L*rho/g)*np.array([[1.0, 0.0],
                                     [0.0, 1.0]])
    #
    cG   = calpha*kG + cbeta*mG
    #
    svar, statev, stress, strann = svarshandl_spring(1, svar, statev, stress, strann)
                
    return  kG , mG , cG, svar, ilf


#####################        NON-LINEAR 1D ROTATIONAL SPRING        ######################

def NL_1DRotspg(iele_disp, coord, par, svar):
    """ Non-Linear 1D Rotational spring element 
    
    Bases on Algorithm of 1-D Rate independent plasticity, Isotropic hardening.
    Chapter 1. Box 1.4, Simo & Hughes. Computational Inelasticity.
    
    Parameters
    ----------
    iele_disp : ndarray
                Nodal global displacement values at time t 
    coord     : ndarray
                Coordinates for the nodes of the element (2, 2).
    par       : ndarray
                Material properties
    svar      : 1D array 
                Element's plastic variables of the increment    
    Returns
    -------
    kG   : ndarray
           Stiffness matrix in global coordinates for the element (2, 2).
    mG   : ndarray
           Mass matrix in global coordinates for the element      (2, 2).
    cG   : ndarray
           Damping matrix in global coordinates for the element   (2, 2).
    svar : 1D array 
           Element's plastic variables of the increment.
    """
    
    Iz     = par[1]
    Emod   = par[2]
    calpha = par[4]
    cbeta  = par[5]
    vec    = coord[1, :] - coord[0, :]
    L      = 1.0 #np.linalg.norm(vec)
    Q      = Q_1Dspring(vec)
    #
    strann = 0.0
    stress = 0.0
    statev = np.zeros([5])
    #
    svar, statev, stress, strann = svarshandl_spring(0, svar, statev, stress, strann)
    strann = strann + (iele_disp[1] - iele_disp[0])*(4.0/L)
    stress, statev, Emod = umat.umat_spring(stress, strann, statev, par)
    #
    #
    kl   = (4*Emod*Iz/L)*np.array([[ 1.0,-1.0],
                                   [-1.0, 1.0]])
    #
    ilf  = np.dot(kl, np.dot(Q,iele_disp))
    kG   = np.dot(np.dot(Q.T, kl), Q)
    #
    mG   = np.zeros_like(kG)
    #
    cG   = calpha*kG + cbeta*mG
    #
    svar, statev, stress, strann = svarshandl_spring(1, svar, statev, stress, strann)
    #      
    return  kG , mG , cG, svar, ilf


#####################        NON-LINEAR 1D PILE P-Y SPRING FOR ROCK        ######################

def svarshandl_1DPileRCK(imode, svars, statev, P, Y):
    """Handles state variables array svars according to imode
    
    Parameters
    ----------
    imode: int
           Handling mode:
           (0) Retrieves state variables
           (1) Updates state variables
    svars: ndarray
           State variables array at all the integration points to be
           retrieved or updated according to imode
    P:     float
           Load at current time/deformation increment.
    Y:     float
           Displacement at current time/deformation increment.
    igp:   int
           Pointer to the svars block for the current integration point.
    Returns
    -------
    svars , statev, stress, strann
    """
    
    if imode == 0:
        P         = svars[0]
        Y         = svars[1]
        statev[0] = svars[2]
        statev[1] = svars[3]

    else:
        svars[0]  = P
        svars[1]  = Y
        svars[2]  = statev[0]
        svars[3]  = statev[1]

     # End if
    return svars, statev, P, Y 

def NL_1DPileRCK(iele_disp, coord, par, svar):
    """ Non-Linear 1D spring element for piles
    
    Contistitive model of a P-y curve for stiff soil.
    Handbook of design of piles under later load. Section 2.9
    
    Parameters
    ----------
    iele_disp : ndarray
                Nodal global displacement values at time t 
    coord     : ndarray
                Coordinates for the nodes of the element (2, 2).
    par       : ndarray
                Material properties
    svar      : 1D array 
                Element's plastic variables of the increment    
    Returns
    -------
    kG   : ndarray
           Stiffness matrix in global coordinates for the element (2, 2).
    mG   : ndarray
           Mass matrix in global coordinates for the element      (2, 2).
    cG   : ndarray
           Damping matrix in global coordinates for the element   (2, 2).
    svar : 1D array 
           Element's plastic variables of the increment.
    """
    vec = coord[1, :] - coord[0, :]
    Q   = Q_1Dspring(vec)

    P = 0.0
    Y = 0.0
    statev = np.zeros([2])
    
    svar, statev, P, Y = svarshandl_1DPileRCK(0, svar, statev, P, Y)
    Y = Y + (iele_disp[1] - iele_disp[0])
    P, statev, E = umat.umat_PileRCK(P, abs(Y), statev, par)
    #
    kl   = E*np.array([[ 1.0,-1.0],
                       [-1.0, 1.0]])
    #
    ilf  = np.dot(kl, np.dot(Q,iele_disp))
    kG   = np.dot(np.dot(Q.T, kl), Q)
    #
    mG   = np.zeros_like(kG)
    cG   = np.zeros_like(kG)
    #
    if Y < 0.0:
       P = P*-1
    #End if   
    svar, statev, P, Y = svarshandl_1DPileRCK(1, svar, statev, P, Y)
            
    return  kG , mG , cG, svar, ilf


##################################################################################################
###################                 LINEAR 1D FRICTIONAL SPRING              #####################

def L_1Dspring(iele_disp, coord, par, svar):
    """ Linear 1D Friccional spring element 
    
    Parameters
    ----------
    iele_disp : ndarray
                Nodal global displacement values at time t 
    coord     : ndarray
                Coordinates for the nodes of the element (2, 2).
    par       : ndarray
                Material properties
    svar      : 1D array 
                Element's plastic variables of the increment    
    Returns
    -------
    kG   : ndarray
           Stiffness matrix in global coordinates for the element (2, 2).
    mG   : ndarray
           Mass matrix in global coordinates for the element      (2, 2).
    cG   : ndarray
           Damping matrix in global coordinates for the element   (2, 2).
    svar : 1D array 
           Element's plastic variables of the increment.
    """
    A      = par[0]
    Emod   = par[2]
    rho    = par[3]
    calpha = par[4]
    cbeta  = par[5]
    g      = par[6]
    vec    = coord[1, :] - coord[0, :]
    L      = np.linalg.norm(vec)
    Q      = Q_1Dspring(vec)
    #
    #
    kl   = (A*Emod/L)*np.array([[ 1.0,-1.0],
                                [-1.0, 1.0]])

    ilf  = np.dot(kl, np.dot(Q,iele_disp)) 
    kG   = np.dot(np.dot(Q.T, kl), Q)
    #
    mG   = (0.5*A*L*rho/g)*np.array([[1.0, 0.0],
                                     [0.0, 1.0]])
    #   
    cG   = calpha*kG + cbeta*mG
    #
    return kG , mG , cG, svar, ilf

##################################################################################################
#################                  LINEAR 1D ROTATIONAL SPRING                  ##################


def L_1DRotspg(iele_disp, coord, par, svar):
    """ Linear 1D Rotational spring element 
    
    Parameters
    ----------
    iele_disp : ndarray
                Nodal global displacement values at time t 
    coord     : ndarray
                Coordinates for the nodes of the element (2, 2).
    par       : ndarray
                Material properties
    svar      : 1D array 
                Element's plastic variables of the increment    
    Returns
    -------
    kG   : ndarray
           Stiffness matrix in global coordinates for the element (2, 2).
    mG   : ndarray
           Mass matrix in global coordinates for the element      (2, 2).
    cG   : ndarray
           Damping matrix in global coordinates for the element   (2, 2).
    svar : 1D array 
           Element's plastic variables of the increment.
    """
    Iz     = par[1]
    Emod   = par[2]
    calpha = par[4]
    cbeta  = par[5]
    vec    = coord[1, :] - coord[0, :]
    L      = np.linalg.norm(vec)
    Q      = Q_1Dspring(vec)
    #
    #
    kl   = (4*Emod*Iz/L)*np.array([[ 1.0,-1.0],
                                   [-1.0, 1.0]])

    ilf  = np.dot(kl, np.dot(Q,iele_disp))
    kG   = np.dot(np.dot(Q.T, kl), Q)
    #
    mG   = np.zeros_like(kG)
    #   
    cG   = calpha*kG + cbeta*mG
    #
    return kG , mG , cG, svar, ilf


##################################################################################################
#################              LINEAR 2D SHEAR - ROTATIONAL SPRING              ##################

def L_2DRotspg(iele_disp, coord, par, svar):
    """ Linear 2D Shear-Rotational spring element 
    
    Parameters
    ----------
    iele_disp : ndarray
                Nodal global displacement values at time t 
    coord     : ndarray
                Coordinates for the nodes of the element (2, 2).
    par       : ndarray
                Material properties
    svar      : 1D array 
                Element's plastic variables of the increment    
    Returns
    -------
    kG   : ndarray
           Stiffness matrix in global coordinates for the element (2, 2).
    mG   : ndarray
           Mass matrix in global coordinates for the element      (2, 2).
    cG   : ndarray
           Damping matrix in global coordinates for the element   (2, 2).
    svar : 1D array 
           Element's plastic variables of the increment.
    """
    Iz     = par[1]
    Emod   = par[2]
    calpha = par[4]
    cbeta  = par[5]    
    vec    = coord[1, :] - coord[0, :]
    L      = np.linalg.norm(vec) 
    Q      = Q_2Dtruss(vec)
    #
    kl = (Iz*Emod/(L*L*L)) * np.array([[ 12.0, 6*L  , -12.0, 6*L  ],
                                       [ 6*L , 4*L*L, -6*L , 2*L*L],
                                       [-12.0, -6*L , 12.0 , -6*L ],
                                       [ 6*L , 2*L*L, -6*L , 4*L*L]])

    ilf  = np.dot(kl, np.dot(Q,iele_disp))
    kG   = np.dot(np.dot(Q.T, kl), Q)
    #
    mG   = np.zeros_like(kG)
    #   
    cG   = calpha*kG + cbeta*mG
    #
    return kG , mG , cG, svar, ilf


##################################################################################################
###################                       LINEAR 2D TRUSS                    #####################

def L_2Dtruss(iele_disp, coord, par, svar):
    """ Linear 2D-2-noded truss element

    Parameters
    ----------
    iele_disp : ndarray
                Nodal global displacement values at time t  
    coord     : ndarray
                Coordinates for the nodes of the element (2, 2).
    par       : ndarray
                Material properties = [Area, Inertia, Emod, Density, Alpha-Damp, Betha-Damp]
    svar      : 1D array 
                Element's plastic variables of the increment      
    Returns
    -------
    kG   : ndarray
           Stiffness matrix in global coordinates for the element (4, 4).
    mG   : ndarray
           Mass matrix in global coordinates for the element      (4, 4).
    cG   : ndarray
           Damping matrix in global coordinates for the element   (4, 4).
    svar : 1D array 
           Element's plastic variables of the increment.
    """
    A      = par[0]
    Emod   = par[2]
    rho    = par[3]
    g      = par[6]
    calpha = par[4]
    cbeta  = par[5]
    vec    = coord[1, :] - coord[0, :]
    L      = np.linalg.norm(vec) 
    Q      = Q_2Dtruss(vec)
    #
    kl = (A*Emod/L) * np.array([[ 1, 0,-1, 0],
                                [ 0, 0, 0, 0],
                                [-1, 0, 1, 0],
                                [ 0, 0, 0, 0]])
    # 
    ilf  = np.dot(kl, np.dot(Q,iele_disp))
    kG = np.dot(np.dot(Q.T, kl), Q)
    #
    ml = (0.5*A*L*rho/g) * np.array([[1.0, 0.0, 0.0, 0.0],
                                     [0.0, 1.0, 0.0, 0.0],
                                     [0.0, 0.0, 1.0, 0.0],
                                     [0.0, 0.0, 0.0, 1.0]])
    mG = ml
    #   
    cG = calpha*kG + cbeta*mG
    #
    return kG , mG , cG, svar, ilf

##################################################################################################
###################                    LINEAR 2D SIMPLE FRAME               ######################
    
def L_sframe(iele_disp, coord, par, svar):
    """ Linear 2D-2-noded frame element
        without axial deformation

    Parameters
    ----------
    iele_disp : ndarray
                Nodal global displacement values at time t
    coord     : ndarray
                Coordinates for the nodes of the element (2, 2).
    par       : ndarray
                Material properties = [Area, Inertia, Emod, Density, Alpha-Damp, Betha-Damp]
    svar      : 1D array 
                Element's plastic variables of the increment      
    Returns
    -------
    kG   : ndarray
           Stiffness matrix in global coordinates for the element (6, 6).
    mG   : ndarray
           Mass matrix in global coordinates for the element      (6, 6).
    cG   : ndarray
           Damping matrix in global coordinates for the element   (6, 6).
    svar : 1D array 
           Element's plastic variables of the increment.
    """
    A      = par[0]
    Iz     = par[1]
    Emod   = par[2]
    rho    = par[3]
    g      = par[6]
    calpha = par[4]
    cbeta  = par[5]    
    vec    = coord[1, :] - coord[0, :]
    L      = np.linalg.norm(vec) 
    Q      = Q_2Dtruss(vec)
    #
    kl = (Iz*Emod/(L*L*L)) * np.array([[ 12.0, 6*L  , -12.0, 6*L  ],
                                       [ 6*L , 4*L*L, -6*L , 2*L*L],
                                       [-12.0, -6*L , 12.0 , -6*L ],
                                       [ 6*L , 2*L*L, -6*L , 4*L*L]])
    #
    ilf  = np.dot(kl, np.dot(Q,iele_disp))
    kG = np.dot(np.dot(Q.T, kl), Q)
    #
    ml = (0.5*A*L*rho/g) * np.array([[1.0, 0.0, 0.0, 0.0],
                                     [0.0, 1.0, 0.0, 0.0],
                                     [0.0, 0.0, 1.0, 0.0],
                                     [0.0, 0.0, 0.0, 1.0]])
    mG = ml
    #
    cG = calpha*kG + cbeta*mG
    #
    return kG , mG , cG, svar, ilf

##################################################################################################
###################                      LINEAR 2D FULL FRAME               ######################
    
def L_fframe(iele_disp, coord, par, svar):
    """ Linear 2D-2-noded complete frame element.
    
    Parameters
    ----------
    iele_disp : ndarray
                Nodal global displacement values at time t
    coord     : ndarray
                Coordinates for the nodes of the element (2, 2).
    par       : ndarray
                Material properties = [Area, Inertia, Emod, Density, Alpha-Damp, Betha-Damp]
    svar      : 1D array 
                Element's plastic variables of the increment      
    Returns
    -------
    kG   : ndarray
           Stiffness matrix in global coordinates for the element (6, 6).
    mG   : ndarray
           Mass matrix in global coordinates for the element      (6, 6).
    svar : 1D array 
           Element's plastic variables of the increment.
    """
    A      = par[0]
    Iz     = par[1]
    Emod   = par[2]
    rho    = par[3]
    g      = par[6]
    calpha = par[4]
    cbeta  = par[5]
    vec    = coord[1, :] - coord[0, :]
    L      = np.linalg.norm(vec)
    Q      = Q_2Dframe(vec)
    #
    kl = (Iz*Emod/(L*L*L)) * np.array([[A*L*L/Iz   , 0.0  , 0.0   , -A*L*L/Iz , 0.0   , 0.0  ],
                                       [ 0.0       , 12.0 , 6*L   , 0.0       , -12.0 , 6*L  ],
                                       [ 0.0       , 6*L  , 4*L*L , 0.0       , -6*L  , 2*L*L],
                                       [-A*L*L/Iz  , 0.0  , 0.0   , A*L*L/Iz  , 0.0   , 0.0  ],
                                       [0.0        ,-12.0 , -6*L  , 0.0       , 12.0  , -6*L ],
                                       [0.0        , 6*L  , 2*L*L , 0.0       , -6*L  , 4*L*L]])
    #
    ilf  = np.dot(kl, np.dot(Q,iele_disp))
    kG = np.dot(np.dot(Q.T, kl), Q)
    #
    ml = (0.5*A*L*rho/g) * np.array([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                     [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                                     [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                                     [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])
    mG =ml
    #
    cG = calpha*kG + cbeta*mG
    #
    return kG , mG , cG, svar, ilf

########################################################################################################
################     LINEAR 3D FULL FRAME BENDING AND SHEAR DEFORMATION EFFECTS      ###################

def L_3Dfframe(iele_disp, coord, par, svar):
    """ Linear 3D-2-noded complete frame element considering bending and shear deformations effects.
    
    Parameters
    ----------
    iele_disp : ndarray
                Nodal global displacement values at time t
    coord     : ndarray
                Coordinates for the nodes of the element (2, 2).
    par       : ndarray
                Material properties = [Area, Inertia, Emod, Density, Alpha-Damp, Betha-Damp]
    svar      : 1D array 
                Element's plastic variables of the increment      
    Returns
    -------
    kG   : ndarray
           Stiffness matrix in global coordinates for the element (6, 6).
    mG   : ndarray
           Mass matrix in global coordinates for the element      (6, 6).
    svar : 1D array 
           Element's plastic variables of the increment.
    """
    A      = par[0]  # Cross sectional area
    J      = par[1]  # Torsional constant
    I2     = par[2]  # Moment of inertia along local axis 2
    I3     = par[3]  # Moment of inertia along local axis 3
    ff2    = par[4]  # Shape factor for effective shear cross area along local axis 2
    ff3    = par[5]  # Shape factor for effective shear cross area along local axis 3
    Emod   = par[6]  # Young's modulus
    Nu     = par[7]  # Poisson's ratio
    rho    = par[8]  # Material's specific weight
    g      = par[9]  # Gravity value
    Rloc   = par[10] # Rotation Angle for local axis 2-3 in degrees
    calpha = par[11] # Rayleight's Damping constant
    cbeta  = par[12] # Rayleight's Damping constant
    #
    vec    = coord[1, :] - coord[0, :]
    L      = np.linalg.norm(vec)
    Q      = Q_3Dframe(vec,Rloc)
    #
    G    = Emod/(2*(1+Nu))
    As2  = A/ff2
    As3  = A/ff3
    Phi2 = 12*Emod*I3/(G*As2*L**2)
    Phi3 = 12*Emod*I2/(G*As3*L**2)
    #
    Kti1 = A*Emod/L
    Kti2 = 12*Emod*I3/((L**3)*(1 + Phi2))
    Kti3 = 12*Emod*I2/((L**3)*(1 + Phi3))
    Kri1 = G*J/L
    Kri2 = 6*Emod*I3/((L**2)*(1 + Phi2))
    Kri3 = 6*Emod*I2/((L**2)*(1 + Phi3))
    Kri4 = (4 + Phi2)*Emod*I3/(L*(1 + Phi2))
    Kri5 = (4 + Phi3)*Emod*I2/(L*(1 + Phi3))
    #
    Krj4 = (2 - Phi2)*Emod*I3/(L*(1 + Phi2))
    Krj5 = (2 - Phi3)*Emod*I2/(L*(1 + Phi3))
    #
    kl = np.array([[    Kti1,     0.0,     0.0,     0.0,     0.0,     0.0, -1*Kti1,     0.0,     0.0,     0.0,     0.0,     0.0],
                   [     0.0,    Kti2,     0.0,     0.0,     0.0,    Kri2,     0.0, -1*Kti2,     0.0,     0.0,     0.0,    Kri2],
                   [     0.0,     0.0,    Kti3,     0.0, -1*Kri3,     0.0,     0.0,     0.0, -1*Kti3,     0.0, -1*Kri3,     0.0],
                   [     0.0,     0.0,     0.0,    Kri1,     0.0,     0.0,     0.0,     0.0,     0.0, -1*Kri1,     0.0,     0.0],
                   [     0.0,     0.0, -1*Kri3,     0.0,    Kri5,     0.0,     0.0,     0.0,    Kri3,     0.0,    Krj5,     0.0],
                   [     0.0,    Kri2,     0.0,     0.0,     0.0,    Kri4,     0.0, -1*Kri2,     0.0,     0.0,     0.0,    Krj4],                  
                   [ -1*Kti1,     0.0,     0.0,     0.0,     0.0,     0.0,    Kti1,     0.0,     0.0,     0.0,     0.0,     0.0],
                   [     0.0, -1*Kti2,     0.0,     0.0,     0.0, -1*Kri2,     0.0,    Kti2,     0.0,     0.0,     0.0, -1*Kri2],
                   [     0.0,     0.0, -1*Kti3,     0.0,    Kri3,     0.0,     0.0,     0.0,    Kti3,     0.0,    Kri3,     0.0],
                   [     0.0,     0.0,     0.0, -1*Kri1,     0.0,     0.0,     0.0,     0.0,     0.0,    Kri1,     0.0,     0.0],
                   [     0.0,     0.0, -1*Kri3,     0.0,    Krj5,     0.0,     0.0,     0.0,    Kri3,     0.0,    Kri5,     0.0],
                   [     0.0,    Kri2,     0.0,     0.0,     0.0,    Krj4,     0.0, -1*Kri2,     0.0,     0.0,     0.0,    Kri4]])
                   
    #
    ilf  = np.dot(kl, np.dot(Q,iele_disp))
    kG   = np.dot(np.dot(Q.T, kl), Q) 
    #
    ml = (0.5*A*L*rho/g) * np.array([[ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                     [ 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                     [ 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                     [ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                     [ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                     [ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                                     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                                     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                                     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                                     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])
    mG =ml
    #
    cG = calpha*kG + cbeta*mG
    #
    return kG , mG , cG, svar, ilf

###################################################################################################

def Q_1Dspring(vec):
    """ 1D-spring trasnformation matrix.
    
    Parameters
    ---------- 
    vec: Vector formed by initial and final node coordenates.
    
    """
    # 
    L  = np.linalg.norm(vec) 
    nx = vec[0]/L
    ny = vec[1]/L
    #
    Q  = np.array([[ nx, ny],
                   [-ny, nx]])
    return Q

def Q_2Dtruss(vec):
    """ 2D-truss trasnformation matrix.
    
    Parameters
    ---------- 
    vec: Vector formed by initial and final node coordenates.
    
    """
    # 
    L  = np.linalg.norm(vec) 
    nx = vec[0]/L
    ny = vec[1]/L
    #
    Q  = np.array([[ nx, ny ,  0,  0],
                   [-ny, nx ,  0,  0],
                   [  0,   0, nx, ny],
                   [  0,   0,-ny, nx]])
    return Q

def Q_2Dframe(vec):
    """ 2D - frame transformation matrix.
    
    Parameters
    ---------- 
    vec: Vector formed by initial and final node coordenates.
    
    """
    # 
    L  = np.linalg.norm(vec) 
    nx = vec[0]/L
    ny = vec[1]/L
    #
    Q  = np.array([[ nx ,  ny ,   0     ,  0   ,  0 , 0 ],
                   [-ny ,  nx ,   0     ,  0   ,  0 , 0 ],
                   [  0 ,  0  ,   1     ,  0   ,  0 , 0 ],
                   [  0 ,  0  ,   0     , nx   , ny , 0 ],
                   [  0 ,  0  ,   0     ,-ny   , nx , 0 ],
                   [  0 ,  0  ,   0     ,  0   , 0  , 1 ]])
    return Q

def Q_3Dframe(vec, Rloc):
    """ 3D - frame transformation matrix.
    
    Parameters
    ---------- 
    vec: Vector formed by initial and final node coordenates.
    
    """
    #
    C1 = np.zeros((3,3))
    C2 = np.zeros((3,3))
    C3 = np.zeros((3,3))
    #
    # Initial reference angles
    Theta = np.arctan2(vec[1], vec[0])                      # Angle between axis 1 and the projection of the element on the plane x-y
    Betha = np.arctan2(vec[2],(vec[0]**2 + vec[1]**2)**0.5) # Angle between new axis 1 and the new direction of the element
    Alpha = Rloc*np.pi/180                                  # Rotation Angle for local axis 2-3 in degrees
    #
    # First rotation
    C1[0,0] = np.cos(Theta)
    C1[0,1] = np.sin(Theta)*-1.0
    C1[1,0] = np.sin(Theta)
    C1[1,1] = np.cos(Theta)
    C1[2,2] = 1.0
    #
    # Second rotation
    C2[0,0] = np.cos(Betha)
    C2[0,2] = np.sin(Betha)*-1.0
    C2[2,0] = np.sin(Betha)
    C2[2,2] = np.cos(Betha)
    C2[1,1] = 1.0
    #
    # Third rotation
    C3[0,0] = 1.0
    C3[1,1] = np.cos(Alpha)
    C3[1,2] = np.sin(Alpha)*-1.0
    C3[2,1] = np.sin(Alpha)
    C3[2,2] = np.cos(Alpha)
    #
    C  = np.dot(np.dot(C3,C2),C1)
    #
    c1 = C[0,0]
    c2 = C[0,1]
    c3 = C[0,2]
    c4 = C[1,0]
    c5 = C[1,1]
    c6 = C[1,2]
    c7 = C[2,0]
    c8 = C[2,1]
    c9 = C[2,2]
    #
    Q  = np.array([[c1,c2,c3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [c4,c5,c6, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [c7,c8,c9, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [ 0, 0, 0,c1,c2,c3, 0, 0, 0, 0, 0, 0],
                   [ 0, 0, 0,c4,c5,c6, 0, 0, 0, 0, 0, 0],
                   [ 0, 0, 0,c7,c8,c9, 0, 0, 0, 0, 0, 0],
                   [ 0, 0, 0, 0, 0, 0,c1,c2,c3, 0, 0, 0],
                   [ 0, 0, 0, 0, 0, 0,c4,c5,c6, 0, 0, 0],
                   [ 0, 0, 0, 0, 0, 0,c7,c8,c9, 0, 0, 0],
                   [ 0, 0, 0, 0, 0, 0, 0, 0, 0,c1,c2,c3],
                   [ 0, 0, 0, 0, 0, 0, 0, 0, 0,c4,c5,c6],
                   [ 0, 0, 0, 0, 0, 0, 0, 0, 0,c7,c8,c9]])
          
    return Q

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    