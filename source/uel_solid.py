# -*- coding: utf-8 -*-
"""
Element subroutines
-------------------
Each UEL subroutine computes the local stiffness matrix for a given
finite element.

New elements can be added by including additional subroutines.

"""
from __future__ import absolute_import, division, print_function
import numpy as np

def uel4nquad_PLK(coord, props , svars , du):
    """Computes the stiffness and mass matrix and consistent internal
       loads vector for a bi-linear quad with a nonlinear
       combined isotropic/kinemtic hardening rate independent
       plasticuty model (see Simo, Juan C., andThomas JR Hughes.
       Computational inelasticity. Vol. 7.
       Springer Science & Business Media, 2006).
    Parameters
    ----------
    coord    : ndarray
             Nodal coordinates.
    props    : ndarray
             Material properties for the element
    svars    : ndarray
             State variables array at all integration points
             4 components of the stress tensor
             4 components of the total strain tensor
             4 components of the elastic strain tensor
             4 components of the plastic strain tensor
             4 components of the back stress tensor
             1 equivalent plastic strain
             1 equivalent stress             
    nsvars   : int
             Number of state variables per integration point.
             4 components of the elastic strain tensor
             4 components of the plastic strain tensor
             4 components of the back stress tensor
             1 equivalent plastic strain
             1 equivalent stress
    du       : ndarray
             Nodal (incremental) displacements vector
             
    Returns
    -------
    kl       : ndarray
             Element stiffness matrix
    ml       : ndarray
             Element mass matrix
    rhsl     : ndarray
             Consistent internal loads vector
    svars    : ndarray
             Updated state variables array for the element

    """
    rho = props[6]
    nstatev = 14
    ntens  = 4
    stress = np.zeros([ntens]) 
    strann = np.zeros([ntens]) 
    dstran = np.zeros([ntens])
    statev = np.zeros([nstatev])
    
    kl     = np.zeros([8, 8])
    ml     = np.zeros([8, 8])
    cl     = np.zeros([8, 8])
    rhsl   = np.zeros([8])
    XW, XP = gpoints2x2()
    
    ngpts = 4
    for i in range(0, ngpts):
        ri  = XP[i, 0]
        si  = XP[i, 1]
        alf = XW[i]        
        igp = i*22
        svars , statev , stress , strann = svarshandl(0 , svars , statev , stress , strann , igp)
        
        ddet, B = stdm4NQ(ri, si, coord)
        dstran  = np.dot(B,du)
        strann  = strann + dstran
        C , stress , statev = umat_PCLK(stress , strann , dstran , statev , props , ntens)
        rhsl = rhsl + np.dot(B.T,stress)*alf*ddet                        
        kl = kl + np.dot(np.dot(B.T,C), B)*alf*ddet
        N  = sha4(ri , si )
        ml = ml + rho*np.dot(N.T , N)*alf*ddet
        
        svars , statev , stress , strann = svarshandl(1 , svars , statev , stress , strann , igp)
          
    return kl, ml, cl, svars, rhsl

def sha4(x, y):
    """Shape functions for a 4-noded quad element

    Parameters
    ----------
    x : float
      x coordinate for a point within the element.
    y : float
      y coordinate for a point within the element.

    Returns
    -------
    N : Numpy array
      Array of interpolation functions.

    Examples
    --------
    We can check evaluating at two different points, namely (0, 0) and
    (1, 1). Thus

    >>> N = sha4(0, 0)
    >>> N_ex = np.array([
    ...    [1/4, 0, 1/4, 0, 1/4, 0, 1/4, 0],
    ...    [0, 1/4, 0, 1/4, 0, 1/4, 0, 1/4]])
    >>> np.allclose(N, N_ex)
    True

    and

    >>> N = sha4(1, 1)
    >>> N_ex = np.array([
    ...    [0, 0, 0, 0, 1, 0, 0, 0],
    ...    [0, 0, 0, 0, 0, 1, 0, 0]])
    >>> np.allclose(N, N_ex)
    True

    """
    N = np.zeros((2, 8))
    H = 0.25*np.array(
        [(1 - x)*(1 - y),
         (1 + x)*(1 - y),
         (1 + x)*(1 + y),
         (1 - x)*(1 + y)])
    N[0, ::2] = H
    N[1, 1::2] = H

    return N



def gpoints2x2():
    """Gauss points for a 2 by 2 grid

    Returns
    -------
    xw : ndarray
      Weights for the Gauss-Legendre quadrature.
    xp : ndarray
      Points for the Gauss-Legendre quadrature.

    """
    xw = np.zeros([4])
    xp = np.zeros([4, 2])
    xw[:] = 1.0
    xp[0, 0] = -0.577350269189626
    xp[1, 0] = 0.577350269189626
    xp[2, 0] = -0.577350269189626
    xp[3, 0] = 0.577350269189626

    xp[0, 1] = 0.577350269189626
    xp[1, 1] = 0.577350269189626
    xp[2, 1] = -0.577350269189626
    xp[3, 1] = -0.577350269189626

    return xw, xp


def svarshandl(imode , svars , statev , stress , strann , igp):
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
        stress = svars[igp:igp+4:1]
        strann = svars[igp+4:igp+8:1]
        statev[0:4:1]   = svars[igp+8:igp+12:1]
        statev[4:8:1]   = svars[igp+12:igp+16:1]
        statev[8:12:1]  = svars[igp+16:igp+20:1]
        statev[12:13:1] = svars[igp+20:igp+21:1]
        statev[13:13:1] = svars[igp+21:igp+22:1]
    else:
        svars[igp:igp+4:1]     = stress
        svars[igp+4:igp+8:1]   = strann
        svars[igp+8:igp+12:1]  = statev[0:4:1]
        svars[igp+12:igp+16:1] = statev[4:8:1]
        svars[igp+16:igp+20:1] = statev[8:12:1]
        svars[igp+20:igp+21:1] = statev[12:13:1]
        svars[igp+21:igp+22:1] = statev[13:14:1]
    
    return svars , statev , stress , strann

def stdm4NQ(r, s, coord):
    """Strain-displacement interpolator B for a 4-noded quad element

    Parameters
    ----------
    r : float
      r component in the natural space.
    s : float
      s component in the natural space.
    coord : ndarray
      Coordinates of the nodes of the element (4, 2).

    Returns
    -------
    ddet : float
      Determinant evaluated at `(r, s)`.
    B : ndarray
      Strain-displacement interpolator evaluated at `(r, s)`.

    """
    nn = 4
    B = np.zeros((4, 2*nn))
    dhdx = 0.25*np.array([
            [s - 1, -s + 1, s + 1, -s - 1],
            [r - 1, -r - 1, r + 1, -r + 1]])
    det, jaco_inv = jacoper(dhdx, coord)
    dhdx = np.dot(jaco_inv, dhdx)
    B[0, ::2] = dhdx[0, :]
    B[1, 1::2] = dhdx[1, :]
    B[3, ::2] = dhdx[1, :]
    B[3, 1::2] = dhdx[0, :]
    return det, B

def jacoper(dhdx, coord):
    """
    Compute the Jacobian of the transformation evaluated at
    the Gauss point

    Parameters
    ----------
    dhdx : ndarray
      Derivatives of the interpolation function with respect to the
      natural coordinates.
    coord : ndarray
      Coordinates of the nodes of the element (nn, 2).

    Returns
    -------
    jaco_inv : ndarray (2, 2)
      Jacobian of the transformation evaluated at `(r, s)`.

    """
    jaco = dhdx.dot(coord)
    det = np.linalg.det(jaco)
    if np.isclose(np.abs(det), 0.0):
        msg = "Jacobian close to zero. Check the shape of your elements!"
        raise ValueError(msg)
    jaco_inv = np.linalg.inv(jaco)
    if det < 0.0:
        msg = "Jacobian is negative. Check your elements orientation!"
        raise ValueError(msg)
    return det, jaco_inv



####
def umat_PCLK(stress , strann , dstran , statev , props , ntens):
    """Return mapping integration algorithm for J2-flow theory
       rate independent plasticity under plane strain conditions
       combined isotropic kinematic hardenig (see Simo, Juan C., and
       Thomas JR Hughes. Computational inelasticity. Vol. 7.
       Springer Science & Business Media, 2006).
    
    Parameters
    ----------
    stress     : ndarray
               Stress tensor at the current integration point at the 
               beginning of the increment.
    strann     : ndarray
               Total strains tensor at the current integration point.
    dstran     : ndarray
               Incremental strains at the current integration point.
    statev     : nd array
               Sate variables array at the current integration point
               at the beginning of the increment.
               4 components of the elastic strain tensor
               4 components of the plastic strain tensor
               4 components of the back stress tensor
               1 equivalent plastic strain
               1 equivalent stress
    props      : ndarray
               Material properties for the element
    ntens      :int
               Number of components for the stress and strain tensors.
               
    Returns
    -------
    
    C          : ndarray
               Constitutive tensor
    stress     : ndarray
               Updated stress tensor at the current integration point.
    statev     : nd array
               Sate variables array at the current integration point
               at the end of the increment.    
    """
    
    toler = 1.0e-7
    istop = 1
    
    eelas   = np.zeros([4])
    eplas   = np.zeros([4])
    xback   = np.zeros([4])
    flow    = np.zeros([4])
    eqplas = 0.0
    smises = 0.0

    statev , eelas , eplas , xback , eqplas , smises = stv_handl( 0 , statev , eelas , eplas , xback , eqplas , smises , ntens)
    
    Emod       = props[0]
    enu        = props[1]
    eg2=Emod/(1.0+enu)
    eg = eg2/2.0    
    sig0    = props[2]
    sigsat  = props[3]
    hrdrate = props[4]
    hmod    = props[5]
       
    C = elas_tensor(enu , Emod)   
    stress = stress + np.dot(C , dstran)
    shydro =(stress[0]+stress[1]+stress[2])/3.0
    eelas  = eelas + dstran
    sdev = deviator(stress)
    stsrel = sdev - xback
    smises = vmises(stsrel)
    
    fbar = np.sqrt(2.0/3.0)*smises    
    syiel0 , syieldk , ehardi , ehardk = uhardnlin(sig0 , sigsat , hrdrate , hmod , eqplas)
    syield   = syiel0
    syieldk0 = syieldk
    
    if fbar > (1.0+toler)*syiel0:
        flow = stsrel/fbar       
        gam_par , eqplas , syieldk , istop = local_NR(syield, syieldk0 , ehardi , ehardk , hmod , sig0 , sigsat , hrdrate , eqplas , fbar , eg)
        for k in range(3):
            xback[k]  = xback[k] + np.sqrt(2.0/3.0)*(syieldk-syieldk0)*flow[k]
            eplas[k]  = eplas[k] + gam_par*flow[k]
            eelas[k]  = eelas[k] - eplas[k]
            stress[k] = flow[k]*syield + xback[k] + shydro
        xback[3] = xback[3] + np.sqrt(2.0/3.0)*(syieldk-syieldk0)*flow[3]
        eplas[3]  = eplas[3] + 2.0*gam_par*flow[3]
        eelas[3]  = eelas[3] - eplas[3]
        stress[3] = flow[3]*syield + xback[3]
        eqplas = eqplas+np.sqrt(2.0/3.0)*gam_par
        C = np.zeros((4 , 4))
        C = plas_tensor(gam_par , fbar , flow , ehardk , ehardi , enu , Emod)       
#
# Store updated values of state variables
#
    statev , eelas , eplas , xback , eqplas , smises = stv_handl( 1 , statev , eelas , eplas , xback , eqplas , smises , ntens)
       
    if istop == 0:
        print('local plasticty algorithm did not converged')
        print('After' , iter , 'iterations')
        print('Last value of the consistency parameter' , gam_par)    
            
    return C , stress , statev

def elas_tensor(nu , E):
    """Assembles the elastic constitutive tensor
    for plane strain idealization
    
    Parameters
    ----------
    
    nu   : float
         Poissons ratio
    E    : float
         Youngs modulus
         
    Returns
    -------
    
    C    : ndarray
         Constitutive tensor

    """
    C = np.zeros((4 , 4))
    opnu   = 1.0+nu
    ominnu = 1.0-2.0*nu
    elam = (E*nu)/(opnu*ominnu)
    eg2=E/(1.0+nu)
    eg=eg2/2.0
    
    
    C[0, 0] = eg2+elam
    C[1, 1] = eg2+elam
    C[2, 2] = eg2+elam   
    C[3, 3] = eg
    
    C[0, 1] = elam
    C[0, 2] = elam   
    C[1, 0] = elam
    C[1, 2] = elam    
    C[2, 0] = elam
    C[2, 1] = elam
    
    return C


def plas_tensor(gam_par , fbar , flow , ehardk , ehardi , nu , E):
    """Assembles the elastoplastic constitutive tensor
    consistent with a return mapping integration algorithm (see Simo and Hughes)
    
    Paramters
    ---------
    
    gam_par    : float
               Consistency parameter
    fbar       : float
               Equivalent stress term in the yield function
    flow       : float
               Yield stress
    ehardk     : float
               Isotropic hardening modulus
    ehardi     : float
               Kinematic hardening modulus
    nu         : float
               Poissons ratio
    E          : float
               Youngs modulus
               
    Return
    ------

    C         : ndarray
              Elastoplastic constitutive tensor           

    """
    C = np.zeros((4 , 4))     
    ebulk3=E/(1.0-2.0*nu)
    eg2=E/(1.0+nu)
    eg=eg2/2.0
    elam=(ebulk3-eg2)/3.0
    eg3 = 3.0*eg
    
    beta1 = 1.0 - eg2*gam_par/fbar
    c1 = 1.0 + (ehardk+ehardi)/eg3
    
    betabar = (1.0/c1)- (eg2*gam_par/fbar)
    effg   = eg*beta1
    effg2  = 2.0*effg
    efflam = elam + (eg2/3.0)*(1-beta1)
    effhrd = -eg2*betabar
    
    for k1 in range(3):
        for k2 in range(3):
            C[k2 , k1] = efflam
        C[k1 , k1] = effg2 + efflam
    C[3 , 3] = effg
    
    for k1 in range(4):
        for k2 in range(4):
            C[k2 , k1] = C[k2 , k1] + effhrd*flow[k2]*flow[k1]
    
    
    return C

def stv_handl(imode , statev , eelas , eplas , xback , eqplas , smises , ntens):
    """Handles state variables according to imode
    Parameters
    ----------
    imode  : int
           Handling mode:
           (0) Retrieves state variables
           (1) Updates state variables
    statev : ndarray
           State variables array.
    eelas  : ndarray
           Elastic strains tensor
    eplas  : ndarray
           Plastic strains tensor
    xback  : ndarray
           Back stress tensor
    eqplas : float
           Equivalent plastic strain
    smises : float
           Equivalent stress
           
    Returns
    -------
    statev , eelas , eplas , xback , eqplas , smises
    according to imode.
    
    """
    
    if imode == 0:
        for k in range(ntens):
            eelas[k]  = statev[k]
            eplas[k]  = statev[k+ntens]
            xback[k]  = statev[k+2*ntens]
        eqplas = statev[3*ntens]
        smises = statev[3*ntens+1]
    else:
        for k in range(ntens):
            statev[k]         = eelas[k]
            statev[k+ntens]   = eplas[k]
            statev[k+2*ntens] = xback[k]
        statev[3*ntens] = eqplas
        statev[3*ntens+1] = smises
    
    return statev , eelas , eplas , xback , eqplas , smises

def deviator(stress):
    """Computes the deviatoric component of
    the stress tensor
    
    Paraemters
    ----------
    stress   : ndarray
             Stress tensor
    
    Returns
    -------
    sdev    : ndarray
            Stress deviator
    
    """
    
    P = proyector()
    sdev = np.dot(P , stress)
    
    return sdev


def local_NR(syield, syieldk0 , ehardi , ehardk , hmod , sig0 , sigsat , hrdrate , eqplas , fbar , eg):
    """Uses a Newton-Raphson approach to find
    the consistency parameter gama_par
    
    Parameters
    ----------
    
    syield    : float
              Yield stress at the current value of alpha_n
    syieldk0  :float
              Kinematic hardening parameter at the current value alpha_n
    ehardi    : float
              Isotropic hardening modulus
    ehardk    : float
              Kinemtic hardening modulus
    hmod      : float
              Kinematic hardening parameter
    sig0      : float
              Initial value of the yield stress
    sigsat    : float
              Saturation stress
    hrdrate   : float
              Isotropic hardening paramter
    eqplas    : float
              Equivalent plastic strain at t_n
    fbar      : float
              Equivalent stress to evaluate the yield function
    eg        : float
              Shear modulus
    
    Returns
    -------
    
    gam_par  : float
             Consistency parameter
    eqplas   : float
             Updated value of the equivalent plastic strain at t_n+1
    syieldk  : float
             Kinematic hardening at t_n+1
    istop    : int
             Convergency index
             1: Convergency
             0: Not convergency

    """
    toler = 1.0e-6
    mxiter = 30
    syieldk = syieldk0
    eg2 = 2.0*eg
    eg3 = 3.0*eg
#
    iter =  1
    gam_par = 0.0
    iflag = 0
    istop = 1
    while iflag == 0:
        f_gamma = -syield+fbar-(eg2*gam_par+np.sqrt(2.0/3.0)*(syieldk-syieldk0))
        teta_2 = 1.0 + ((ehardi+ehardk)/eg3)
        fjac = -eg2*teta_2
        gam_par = gam_par - f_gamma/fjac
        eqplas = eqplas + np.sqrt(2.0/3.0)*gam_par
        syield , syieldk , ehardi , ehardk = uhardnlin(sig0 , sigsat , hrdrate , hmod , eqplas)
        if np.abs(f_gamma/fjac) < toler:
            iflag = 1
        else:           
            iter = iter + 1
            if iter>mxiter:
                iflag = 1
                istop = 0
                print('Local N-R diverged')
                print('Last found increment' , np.abs(f_gamma/fjac))
                print('Last found parameter' , gam_par)


    return gam_par , eqplas , syieldk , istop

def uhardnlin(sig0 , sigsat , hrdrate , hmod , eqplas):
    """Evaluates hardening at the specified value of the
       equivalent plastic strain
       
    Parameters
    ----------
    sig0      : float
              Initial value of the yield stress
    sigsat    : float
              Saturation stress
    hrdrate   : float
              Isotropic hardening paramter
    hmod      : float
              Kinematic hardening parameter
    eqplas    : float
              Equivalent plastic strain at t_n
              
    Returns
    -------
    syieldi   : float
              Yield stress at the current value of alpha_n
    syieldk   : float
              Kinematic hardening parameter at the current value alpha_n    
    ehardi    : float
              Isotropic hardening modulus
    ehardk    : float
              Kinemtic hardening modulus    

    """
    syieldi = np.sqrt(2.0/3.0)*(sig0 + sigsat*(1.0 - np.exp(-hrdrate*eqplas)))
    ehardi  = sigsat*hrdrate*(np.exp(-hrdrate*eqplas))
    syieldk = hmod*eqplas
    ehardk = hmod


    return syieldi , syieldk , ehardi , ehardk
    
def proyector():
    """Assembles array P which extracts the
    deviatoric component of a symmetric 2D
    tensor

    """
    P = np.zeros((4 , 4))
    P[0,0] = 2.0/3.0
    P[1,1] = 2.0/3.0
    P[2,2] = 2.0/3.0
    P[3,3] = 1.0
    
    P[0,1] = -1.0/3.0
    P[0,2] = -1.0/3.0
    P[1,0] = -1.0/3.0
    P[1,2] = -1.0/3.0
    P[2,0] = -1.0/3.0
    P[2,1] = -1.0/3.0
    
    return P

def vmises(s):
    """Computes the Von Misses stress for s
    
    Parameters
    ----------
    
    s    : ndarray
         Stress tensor
         
    Returns
    -------
    
    svm : float
        Von Misses stress

    """
    svm = (s[0]-s[1])**2 + (s[1]-s[2])**2 + (s[2]-s[0])**2 + 6.0*s[3]**2
    svm = np.sqrt(svm/2.0)

    return svm 


if __name__ == "__main__":
    import doctest
    doctest.testmod()
