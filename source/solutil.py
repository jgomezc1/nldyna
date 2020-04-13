"""
solutil.py
----------

Utilities for solution of FEM systems

"""
from __future__ import division, print_function
import numpy as np
from scipy.sparse.linalg import splu
import assemutil as ass
import preprocesor as pre
import copy
from numpy import ndarray
from numpy.linalg import solve
from scipy.sparse.csr import csr_matrix
from scipy.sparse.linalg import spsolve
from datetime import datetime

def system_sol(inipar, Up, Vp, neq, RHSG, MG, KG, CG, DME, Msvar, ILF, elements, mats, nodes, ac, const, cst_neq):
    """
    System's nonlinear dynamic or static solution
    
    Parameters:
    ----------
    inipar   = Global parameters
    Up       = Initial nodal displacements vector
    Vp       = Initial nodal velocities vector
    neq      = Number of free DOF
    RHSG     = Global loads assembly
    MG       = Global mass matrix
    KG       = Global stiffness matrix
    CG       = Global damping matrix
    DME      = Global index aseembly operator
    Msvar    = Global state variables for an increment
    ILF      = Internal local forces of the elements
    elements = Global array of elements
    mats     = Global array of materials profiles
    nodes    = Global array of nodes
    ac       = Integration method's constants (ONLY for NLDYNA)
    const    = Nodal constraints information
    cst_neq  = Number of equations after constraints
    
    Return:
    -------
    U        = System's displacement response
    MvarsGen = History of state variables of the system solution      
    
    """
    NLSTA  = int(inipar[5,0])
    NLDYNA = int(inipar[6,0])
    #
    ninc, T, dt, ac, theta = pre.intparams(inipar)
    Tol                    = inipar[2,0]
    maxite                 = inipar[3,0]
    #
    if (NLSTA == 1) and (NLDYNA == 0):
         U = np.zeros((cst_neq, ninc))
         start_time          = datetime.now()
         U, MvarsGen, ILFGen = NonLinResponseST(inipar, DME, Msvar, ILF, elements, mats, nodes, neq, ninc, dt, U, RHSG, KG, Tol, maxite, ac, const)
         end_time            = datetime.now()
         print('Duration for system solution: {}'.format(end_time - start_time))
         #
    elif (NLSTA == 0) and (NLDYNA == 1):
         KE      = ass.effective(KG, MG, CG, ac)
         U, V, A = inicond(Up, Vp, ninc , cst_neq, RHSG, MG, KG, CG)
         print("Finished initial conditions....: {}".format(0))
         del KG
         start_time          = datetime.now()
         U, MvarsGen, ILFGen = NonLinResponse(inipar, DME, Msvar, ILF, elements, mats, nodes, neq, ninc, dt, theta, ac, U, V, A, RHSG, MG, CG, KE, Tol, maxite, const)
         end_time            = datetime.now()
         print('Duration for system solution: {}'.format(end_time - start_time))
         #
    elif (NLSTA == 1) and (NLDYNA == 1):
         inipar[5,0] = 0
         KE      = ass.effective(KG, MG, CG, ac)
         U, V, A = inicond(Up, Vp, ninc , cst_neq, RHSG, MG, KG, CG)
         print("Finished initial conditions....: {}".format(0))
         del KG
         start_time          = datetime.now()
         U, MvarsGen, ILFGen = NonLinResponse(inipar, DME, Msvar, ILF, elements, mats, nodes, neq, ninc, dt, theta, ac, U, V, A, RHSG, MG, CG, KE, Tol, maxite, const)
         end_time            = datetime.now()
         print('Duration for system solution: {}'.format(end_time - start_time))
         #
    # End if 
    return U, MvarsGen, ILFGen 

def static_sol(mat, rhs):
    """Solve a static problem [mat]{u_sol} = {rhs}
    """
    if type(mat) is csr_matrix:
        u_sol = spsolve(mat, rhs)
    elif type(mat) is ndarray:
        u_sol = solve(mat, rhs)
    else:
        raise Exception("Not supported matrix storage scheme!")

    return u_sol

def inicond_U_V(neq):
    """
     Currently homogeneous initial displacement and velocity conditions only
     
     Parameters:
     -----------
     - neq  = number of equations (DOF)

    """
    Up  =np.zeros([neq],dtype=np.float)
    Vp  =np.zeros([neq],dtype=np.float)
    
    return Up, Vp

def inicond(Up, Vp, ninc , neq, RHSG , MG , KG , CG ):
    """
     Currently homogeneous initial conditions only
     
     Parameters:
     -----------
     - Up = Initial nodal displacements
     - Vp = Initial nodal velocities 
     - ninc = number of time increments
     - neq  = number of equations (DOF)
     - RHSG = Nodal time dependent loads
     - MG = Global mass matrix
     - KG = Global stiffness matrix
     - CG = Global damping matrix
     
    """
    Ap  =np.zeros([neq  ],dtype=np.float)
    U   =np.zeros([neq,ninc],dtype=np.float)
    V   =np.zeros([neq,2],dtype=np.float)
    A   =np.zeros([neq,2],dtype=np.float)
    F   =np.zeros([neq],dtype=np.float)
    FE  =np.zeros([neq],dtype=np.float)
    F = RHSG[:, 0]
    FS = KG.dot(Up)
    FD = CG.dot(Vp)
    FE = F - FD - FS
    #
    Ap = static_sol(MG,FE)   
    A[:, 0] = Ap
    #
    return U , V , A

def NonLinResponseST(inipar, DME, Msvar, ILF, elements, mats, nodes, neq, ninc, dt, U, F, KG, Tol, maxite, ac, const):
    """
    Uses the incremental formulation of Wilson theta method to perform implicit time integration (Chopra, Table 15.3.3)
    The method uses Gamma = 1/2 and Betha = 1/6 to perform the time integration
    
    Parameters
    ----------
    DME      : Assembly operator that indicates degrees of freedom correspondence.
    Msvar    : List of list of the state variables of each element at the last converged increment.
    elements : Array with the number for the nodes in each element, and elastoplasticity type.
    mats     : Array with the material/section definitions.
    nodes    : Array with the nodal numbers and coordinates.
    neq      : Integer. Number of equations (degrees of freedom)
    ninc     : Integer. Number of time increments
    dt       : Float. Time step
    U        : Float arrays with the nodal point displacement
    F        : Float array with the point load amplitudes.
    KG       : Float array. Stiffness matrix.
    Tol      : Error tolerance for the newton-raphson iterations.
    maxite   : Maximum number of iterations allowed
    ac       : Integration Method constants (Not used)
    const    : Nodal constraint information
    
    """
    #
    MvarsGen = []
    ILFGen   = []
    Msvar_flag = inipar[14]
    ILF_flag = inipar[15]
    #
    for k in range(ninc-1): 
        dtaFE = np.zeros([neq], dtype=np.float)
        dtaFE = (F[:,k+1] - F[:,k])
        #        
        U, KG, Msvar, ILF = IncrementResponseST(inipar, DME, Msvar, ILF, elements, mats, nodes, neq, KG, dtaFE, U, dt, k, Tol, maxite, ac, const)
        #
        if Msvar_flag == 1:
           Svar = copy.deepcopy(Msvar)
           MvarsGen.append(Svar)
        # End if
        if ILF_flag == 1:
           Iforces = copy.deepcopy(ILF)
           if k != 0:
             for i in range (len(Msvar)):
                 Iforces[i] = Iforces[i] + ILFGen[k-1][i]
             # End for i
           # End if  
           ILFGen.append(Iforces)
        # End if
    #End for
     
    return U, MvarsGen, ILFGen

def IncrementResponseST(inipar, DME, Msvar, ILF, elements, mats, nodes, neq, KG, dFext, U, dtaT, k, Tol, maxite, ac, const):
    """
    This function uses the Newton Raphson scheme to compute displacement, acceleration and velocity
    vectors of elastoplastic systems in a time increment of the MDOF system. This function uses the 
    theta wilson time integration scheme.
    
    Parameters
    ----------
    DME      : Assembly operator that indicates degrees of freedom correspondence.
    Msvar    : List of list of the state variables of each element at the last converged increment.
    elements : Array with the number for the nodes in each element, and elastoplasticity type.
    mats     : Array with the material/section definitions.
    nodes    : Array with the nodal numbers and coordinates. 
    neq      : Number of equations (degrees of freedom)
    KG       : Float array. Stiffness matrix.
    Fext     : Float array with nodal external forces.
    U        : Float array with the nodal point displacements.
    dtaT     : Float. Time step.
    k        : time increment of the integration.
    Tol      : Error tolerance for the newton-raphson iterations
    maxite   : Maximum number of iterations allowed
    ac       : Integration constant (not used)
    const    : Nodal constraint information
    
    """
    LU = splu(KG)
    dUp = LU.solve(dFext)
    # Internal forces
    KG, MG, CG, Msvar, dFint, ILF = ass.assembler(inipar, dUp, Msvar, ILF, elements, mats, nodes, neq, DME, ac, const)
    #
    ResF = (np.linalg.norm(dFext) - np.linalg.norm(dFint))
    dtaR = dFext - dFint
    #
    solFlag = 0
    nite    = 0
    #
    if abs(ResF) < Tol:
       solFlag = 1 # Elastic increment
    else:
         while solFlag == 0:  
           LU   = splu(KG)
           dtaU = LU.solve(dtaR)
           dUp  = dUp + dtaU
           #     
           KG, MG, CG, Msvar, dFint_i, ILF = ass.assembler(inipar, dUp, Msvar, ILF, elements, mats, nodes, neq, DME, ac, const)
           #
           ResF = (np.linalg.norm(dFext) - np.linalg.norm(dFint_i))
           dtaR = dFext - dFint_i
           #
           ResdU = np.dot(dtaR, dtaU) / np.dot(dFext, dUp)
           #
           nite = nite + 1
           if (abs(ResF) < Tol) and (ResdU < Tol):
              solFlag = 1
              print('Convergency reached after ',nite,' iterations', 'at increment ',k,' (',round(dtaT*k,4),'sec)')
           else:
              if nite > maxite:
                 solFlag = 1
                 print('Maximum number of iterations allowed has been reached at increment = ', k)
              # End if
           # End if
         # End while
    # End if 
    # 
    U[:,k+1] = U[:,k] + dUp 
    #
    return U, KG, Msvar, ILF


def NonLinResponse(inipar, DME, Msvar, ILF, elements, mats, nodes, neq, m, dt, theta, ac, U, V, A, F, MG, CG, KE, Tol, maxite, const):
    """
    Uses the incremental formulation of Wilson theta method to perform implicit time integration (Chopra, Table 15.3.3)
    The method uses Gamma = 1/2 and Betha = 1/6 to perform the time integration
    
    Parameters
    ----------
    DME      : Assembly operator that indicates degrees of freedom correspondence.
    Msvar    : List of list of the state variables of each element at the last converged increment.
    elements : Array with the number for the nodes in each element, and elastoplasticity type.
    mats     : Array with the material/section definitions.
    nodes    : Array with the nodal numbers and coordinates.
    neq      : Integer. Number of equations (degrees of freedom)
    m        : Integer. Number of time increments
    dt       : Float. Time step.
    ass      : Float array. Integration constants.
    theta    : Float. Integration parameter.
    U, V, A  : Float arrays with the nodal point displacement,
               velocity and acceleration. Must be passed with the
               initial conditions.
    F        : Float array with the point load amplitudes.
    MG       : Float array. Mass matrix.
    CG       : Float array. Dammping matrix.
    KE       : Float array. Effective stiffness matrix.
    Tol      : Error tolerance for the newton-raphson iterations.
    maxite   : Maximum number of iterations allowed
    const    : Nodal constraints information
    
    """
    #
    MvarsGen = []
    ILFGen = []
    Msvar_flag = inipar[14]
    ILF_flag = inipar[15]
    #
    for k in range(m-1): 
        VV    = np.zeros([neq],dtype=np.float)
        AA    = np.zeros([neq],dtype=np.float)
        dtaFe = np.zeros([neq],dtype=np.float)
        dtaFE = np.zeros([neq],dtype=np.float)
        #
        dtaFe = theta*(F[:,k+1] - F[:,k])
        a     = ac[0]*MG + 3.0*CG
        b     = 3.0*MG + ac[1]*CG 
        #
        VV    = V[:,0]
        AA    = A[:,0]
        dtaFD = a.dot(VV) 
        dtaFI = b.dot(AA)
        dtaFE = dtaFe + dtaFI + dtaFD
        #        
        U, A, V, KE, Msvar, ILF = IncrementResponse(inipar, DME, Msvar, ILF, elements, mats, nodes, neq, KE, dtaFE, A, V, U, ac, dt, k, Tol, maxite, const)
        #
        A[:,0] = A[:,1]
        V[:,0] = V[:,1]
        #
        if Msvar_flag == 1:
           Svar = copy.deepcopy(Msvar)  
           MvarsGen.append(Svar)
        # End if
        if ILF_flag == 1:
           Iforces = copy.deepcopy(ILF)  
           if k != 0:
             for i in range (len(Msvar)):
                 Iforces[i] = Iforces[i] + ILFGen[k-1][i]
             # End for i
           # End if   
           ILFGen.append(Iforces)
        # End if
    #End for
     
    return U, MvarsGen, ILFGen

def IncrementResponse(inipar, DME, Msvar, ILF, elements, mats, nodes, neq, KE, dFext, A, V, U, ac, dtaT, k, Tol, maxite, const):
    """
    This function uses the Newton Raphson scheme to compute displacement, acceleration and velocity
    vectors of elastoplastic systems in a time increment of the MDOF system. This function uses the 
    theta wilson time integration scheme.
    
    Parameters
    ----------
    DME      : Assembly operator that indicates degrees of freedom correspondence.
    Msvar    : List of list of the state variables of each element at the last converged increment.
    elements : Array with the number for the nodes in each element, and elastoplasticity type.
    mats     : Array with the material/section definitions.
    nodes    : Array with the nodal numbers and coordinates. 
    neq      : Number of equations (degrees of freedom)
    KE       : Float array. Effective stiffness matrix.
    Fext     : Float array with nodal external forces.
    A        : Float arrays with the nodal point acceleration.
    V        : Float arrays with the nodal point velocity.
    U        : Float array with the nodal point displacements.
    ac       : Float array. Integration constants.
    dtaT     : Float. Time step.
    k        : time increment of the integration.
    Tol      : Error tolerance for the newton-raphson iterations
    maxite   : Maximum number of iterations allowed
    const    : Nodal constraints information
    
    """
    LU = splu(KE)
    dUp = LU.solve(dFext)
    # Internal forces
    KG, MG, CG, Msvar, dFint, ILF = ass.assembler(inipar, dUp, Msvar, ILF, elements, mats, nodes, neq, DME, ac, const)
    KE = ass.effective(KG, MG, CG, ac)
    #
    ResF = np.linalg.norm(dFext) - np.linalg.norm(dFint)
    dtaR = dFext - dFint
    #
    solFlag = 0
    nite    = 0
    #
    if abs(ResF) < Tol:
       solFlag = 1 # Elastic increment
    else:
         while solFlag == 0:
           LU   = splu(KE)
           dtaU = LU.solve(dtaR)
           dUp  = dUp + dtaU
           #     
           KG, MG, CG, Msvar, dFint_i, ILF = ass.assembler(inipar, dUp, Msvar, ILF, elements, mats, nodes, neq, DME, ac, const)
           KE = ass.effective(KG, MG, CG, ac)
           #
           ResF = (np.linalg.norm(dFext) - np.linalg.norm(dFint_i))
           dtaR = dFext - dFint_i
           #
           ResdU = np.dot(dtaR, dtaU) / np.dot(dFext, dUp)
           #
           nite = nite + 1
           if (abs(ResF) < Tol) and (ResdU < Tol):
              solFlag = 1
              print('Convergency reached after ',nite,' iterations', 'at increment ',k,' (',round(dtaT*k,4),'sec)')
           else:
              if nite > maxite:
                 solFlag = 1
                 print('Maximum number of iterations allowed has been reached at increment = ', k)
              # End if
           # End if
         # End while
    # End if 
    # 
    dAp = ac[3]*dUp - ac[0]*V[:,0] - 3.0*A[:,0]
    dA  = ac[4]*dAp
    dV  = dtaT*A[:,0] + ac[5]*dA
    dU  = dtaT*V[:,0] + ac[6]*A[:,0] + ac[7]*dA
    #
    A[:,1]   = A[:,0] + dA
    V[:,1]   = V[:,0] + dV
    U[:,k+1] = U[:,k] + dU 
    #
    #
    return U, A, V, KE, Msvar, ILF