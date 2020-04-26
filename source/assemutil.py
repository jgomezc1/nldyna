# -*- coding: utf-8 -*-
"""
assemutil.py
------------

Functions to assemble the system of equations for the Finite Element
Analysis.

"""
from __future__ import division, print_function
import numpy as np
import uelutil as ue
import femutil as fem
import uel_solid as slds
import constutil as cst

def eqcounter(nodes):
    """Counts active equations and creates BCs array IBC

    Parameters
    ----------
    nodes : ndarray
      Array with nodes coordinates and boundary conditions.

    Returns
    -------
    neq : int
      Number of equations in the system after removing the nodes
      with imposed displacements.
    IBC : ndarray (int)
      Array that maps the nodes with number of equations.

    """
    nnodes = nodes.shape[0]
    IBC = np.zeros([nnodes, 6], dtype=np.integer)
    IBC[:,:] = nodes[:,4:]
    neq = 0
    for i in range(nnodes):
        for j in range(6):
            if IBC[i, j] == 0:
                IBC[i, j] = neq
                neq = neq + 1
            # End if
        # End for j
    # End for i
    return neq, IBC


def DME(nodes, elements):
    """Counts active equations, creates BCs array IBC[]
    and the assembly operator DME[]

    Parameters
    ----------
    nodes    : ndarray.
      Array with the nodal numbers and coordinates.
    elements : ndarray
      Array with the number for the nodes in each element.
    iet :  int (This is not an input for this function)
      Type of element. These are:
        0.  Linear - 1D Spring in X direction.   
        1.  Linear - Simple frame.
        2.  Linear - Full frame.
        3.  Linear - 2D Truss.
        4.  NonLin - 4 noded plate.
        5.  NonLin - 1D Spring in X direction.
        6.  Linear - 2D Shear-Rotational spring.
        7.  Linear - 1D Rotational spring.
        8.  NonLin - 1D Rotational spring.
        9.  NonLin - 1D pile spring P-Y curve for stiff soil
        10. Linear - 3D full frame bending and shear effects

    Returns
    -------
    DME : ndarray (int)
      Assembly operator.
    IBC : ndarray (int)
      Boundary conditions array.
    neq : int
      Number of active equations in the system.

    """
    nels = len(elements)
    IELCON = np.zeros([nels, 4], dtype=np.integer) # 4 por que es el máximo numero de nodos que conforman un elemento
    DME = []
    
    neq, IBC = eqcounter(nodes)

    for i in range(nels):
        dme = []
        iet = elements[i][1]
        ndof, nnodes, ngpts = fem.eletype(iet)
        for j in range(nnodes):
            IELCON[i, j] = elements[i][j+3]
            kk = IELCON[i, j]
            #
            if iet == 0:
               dme.append(IBC[kk, 0])
            elif iet == 1:
               dme.append(IBC[kk, 0])
               dme.append(IBC[kk, 1])
               dme.append(IBC[kk, 5])
            elif iet == 2:
               dme.append(IBC[kk, 0])
               dme.append(IBC[kk, 1])
               dme.append(IBC[kk, 5])
            elif iet == 3:
               dme.append(IBC[kk, 0])
               dme.append(IBC[kk, 1])
            elif iet == 4:
               dme.append(IBC[kk, 0])
               dme.append(IBC[kk, 1]) 
            elif iet == 5:
               dme.append(IBC[kk, 0])
            elif iet == 6:
               dme.append(IBC[kk, 1])
               dme.append(IBC[kk, 5])
            elif iet == 7:
               dme.append(IBC[kk, 5])
            elif iet == 8:
               dme.append(IBC[kk, 5])
            elif iet == 9:
               dme.append(IBC[kk, 0])
            elif iet == 10:
               dme.append(IBC[kk, 0])
               dme.append(IBC[kk, 1]) 
               dme.append(IBC[kk, 2]) 
               dme.append(IBC[kk, 3]) 
               dme.append(IBC[kk, 4]) 
               dme.append(IBC[kk, 5]) 
            #End if
        # End for j
        DME.append(dme)
    # ENd for i
    return DME, IBC, neq

def eleDisp(Up, iet, dme, i):
    """
    This function find displacements of element i. Returns 1D array
    with x, y and rotation displacements values in each node. 
    
    Parameters
    ----------
    Up     : Displacement al time t + dtaT (not converged).
    iet    : Element i type
    dme    : Assembly operator that indicates degrees of freedom correspondence of the element.
    i      : Element indicator 
    
    """
    ndof, nnodes, ngpts = fem.eletype(iet)
    ele_disp = np.zeros((ndof,))
    #
    for j in range (ndof):
         ID_dof = dme[j]
         if ID_dof == -1:
              ele_disp[j] = 0.0
         else:
              ele_disp[j] = Up[ID_dof]
         # End if
    #End for
    return ele_disp


def retriever(Up, svar, dme, ac, elements , mats , nodes , i, uel=None):
    """Computes the elemental stiffness matrix of element i

    Parameters
    ----------
    Up       : ndarray
               Array with the nodal displacement at time t
    svar     : python list
               List with the state variables of the element of the last converged increment.
    dme      : ndarray (int)
               Assembly operator that indicates degrees of freedom correspondence of the element.
    ac       : 1D array
               Integration constants
    elements : ndarray
               Array with the number for the nodes in each element.
    mats     : ndarray.
               Array with the material profiles.
    nodes    : ndarray.
               Array with the nodal numbers and coordinates.
    i        : int.
               Identifier of the element to be assembled.

    Returns
    -------
    kloc : ndarray (float)
           Array with the local stiffness matrix.
    mloc : ndarray (float)
           Array with the local mass matrix.
    cloc : ndarray (float)
           Array with the local damping matrix.
    svar : Python list
           List with the state variables of the current element.
    ndof : int.
           Number of degrees of freedom of the current element.           
    """
    #
    iet    = elements[i][1]
    #
    ndof, nnodes, ngpts = fem.eletype(iet)
    iele_disp = eleDisp(Up, iet, dme, i)
    #
    elcoor = np.zeros([nnodes, 3])
    im     = np.int(elements[i][2])
    par    = mats[im]
    
    for j in range(nnodes):
        IELCON       = elements[i][j+3]  # Index of each node of the element
        elcoor[j, 0] = nodes[IELCON, 1]  # X coordinate of the node
        elcoor[j, 1] = nodes[IELCON, 2]  # X coordinate of the node
        elcoor[j, 2] = nodes[IELCON, 3]  # X coordinate of the node
    # End for j
    #
    if uel is None:
        if   iet == 0:
             kloc , mloc , cloc , svar, ilf = ue.L_1Dspring(iele_disp, elcoor, par, svar)
        if   iet == 1:
             kloc , mloc , cloc , svar, ilf = ue.L_sframe(iele_disp, elcoor, par, svar)
        elif iet == 2:
             kloc , mloc , cloc , svar, ilf = ue.L_fframe(iele_disp, elcoor, par, svar)
        elif iet == 3:
             kloc , mloc , cloc , svar, ilf = ue.L_2Dtruss(iele_disp, elcoor, par, svar)
        elif iet == 4: 
             kloc , mloc , cloc , svar, ilf = slds.uel4nquad_PLK(elcoor, par, svar, iele_disp)
        elif iet == 5: 
             kloc , mloc , cloc , svar, ilf = ue.NL_1Dspring(iele_disp, elcoor, par, svar)
        elif iet == 6: 
             kloc , mloc , cloc , svar, ilf = ue.L_2DRotspg(iele_disp, elcoor, par, svar)
        elif iet == 7: 
             kloc , mloc , cloc , svar, ilf = ue.L_1DRotspg(iele_disp, elcoor, par, svar)
        elif iet == 8: 
             kloc , mloc , cloc , svar, ilf = ue.NL_1DRotspg(iele_disp, elcoor, par, svar)
        elif iet == 9: 
             kloc , mloc , cloc , svar, ilf = ue.NL_1DPileRCK(iele_disp, elcoor, par, svar) 
        elif iet == 10:
             kloc , mloc , cloc , svar, ilf = ue.L_3Dfframe(iele_disp, elcoor, par, svar)     
    else:
        kloc, ndof, iet = uel(elcoor, par)
    
    return kloc, mloc, cloc, svar, ndof, ilf


def assembler(inipar, Up, Msvar, ILF, elements, mats, nodes, neq, DME, ac, const, sparse=False, uel=None):
    """Assembles the global stiffness matrix

    Parameters
    ----------
    Up       : ndarray (float)
               Nodal displacement values at time t
    Msvar    : Python list
               List of list with the state variables of all the elements at the last converged increment.
    ILF      : Python list
               List of list with the history of internal local forces of each element 
    elements : ndarray (int)
               Array with the number for the nodes in each element.
    mats     : ndarray (float)
               Array with the material profiles.
    nodes    : ndarray (float)
               Array with the nodal numbers and coordinates.
    DME      : ndarray (int)
               Assembly operator.
    neq      : int
               Number of active equations in the system.
    ac       : 1D array
               Integration constants
    const    : ndarray (int)
               Nodal constraints 
    sparse   : boolean (optional)
               Boolean variable to pick sparse assembler. It is True by default.
    uel      : callable function (optional)
               Python function that returns the local stiffness matrix.

    Returns
    -------
    KG : ndarray (float)
      Array with the global stiffness matrix. It might be
      dense or sparse, depending on the value of _sparse_
    MG : ndarray (float)
      Array with the global mass matrix. It might be
      dense or sparse, depending on the value of _sparse_
    CG : ndarray (float)
      Array with the global damping matrix. It might be
      dense or sparse, depending on the value of _sparse_
    IGF : ndarray (float)
      Array with the internal global forces. It might be
      dense or sparse, depending on the value of _sparse_
    ILF : ndarray (float)
      List with the internal local forces.
    """
    #
    if sparse:
        print('Sparse assembly method is not incluyed yet!')
    else:
        KG, MG, CG, Msvar, IGF, ILF = dense_assem(inipar, Up, Msvar, elements, mats, nodes, neq, DME, ac, ILF, const, uel=uel)
     #End if
    
    return KG, MG, CG, Msvar, IGF, ILF

def dense_assem(inipar, Up, Msvar, elements, mats, nodes, neq, DME, ac, ILF, const, uel=None):
    """
    Assembles the global stiffness matrix _KG_
    using a dense storing scheme

    Parameters
    ----------
    Up       : ndarray (float)
               Nodal displacement values at time t
    Msvar    : Python list
               List of list with the state variables of all the elements at the last converged increment.
    elements : ndarray (int)
               Array with the number for the nodes in each element.
    mats     : ndarray (float)
               Array with the material profiles.
    nodes    : ndarray (float)
               Array with the nodal numbers and coordinates.
    DME      : ndarray (int)
               Assembly operator.
    neq      : int
               Number of active equations in the system.
    ac       : 1D array
               Integration constants
    ILF      : Python list
               List of list with the history of internal local forces of each element
    const    : List
               Nodal constraints
    uel      : callable function (optional)
               Python function that returns the local stiffness matrix.

    Returns
    -------
    KG : ndarray (float)
         Array with the global stiffness matrix in a dense numpy array.
    MG : ndarray (float)
      Array with the global mass matrix in a dense numpy array.
    CG : ndarray (float)
      Array with the global damping matrix in a dense numpy array.
      
    """
    
    NLSTA  = int(inipar[5,0])
    NLDYNA = int(inipar[6,0])
    #
    DPH_flag = const[0]
    CST_flag = const[1]
    #
    if (NLSTA == 1) and (NLDYNA == 0):
         DYN_flag = 0
         #
    elif (NLSTA == 0) and (NLDYNA == 1):
         DYN_flag = 1
         #
    elif (NLSTA == 1) and (NLDYNA == 1):
         DYN_flag = 1
         #
    # End if  

    KG  = np.zeros((neq, neq))
    MG  = np.zeros((neq, neq))
    CG  = np.zeros((neq, neq))
    IGF = np.zeros(neq)
    #
    nels = len(elements)
    #
    if (DPH_flag != 0) or (CST_flag != 0):
        if (int(inipar[16,0]) != 0):
            Upn = cst.disp_retrieve(neq, const, Up)
        else:
            Upn = Up 
        #End if
    else:
        Upn = Up
    # End if     
    #
    for el in range(nels):
        kloc, mloc, cloc, svar, ndof, ilf = retriever(Upn, Msvar[el][:], DME[el][:], ac, elements, mats, nodes, el, uel=uel)
        dme          = DME[el][:]
        Msvar[el][:] = svar[:]
        ILF[el][:]   = ilf[:]
        for row in range(ndof):
            glob_row = dme[row]
            if glob_row != -1:
                for col in range(ndof):
                    glob_col = dme[col]
                    if glob_col != -1:
                        KG[glob_row, glob_col] = KG[glob_row, glob_col] +\
                                                 kloc[row, col]                      
                        #
                        if DPH_flag == 1:
                             MG[glob_row, glob_col] = MG[glob_row, glob_col] +\
                                                      mloc[row, col]
                        if DYN_flag == 1:                         
                             MG[glob_row, glob_col] = MG[glob_row, glob_col] +\
                                                      mloc[row, col]
                             CG[glob_row, glob_col] = CG[glob_row, glob_col] +\
                                                      cloc[row, col]
                        #End if 
                    #End if
                 #End for col
              #End if
          #End for row
    #End for el         
    #
    if (DPH_flag != 0) or (CST_flag != 0):
         if (int(inipar[16,0]) == 0):
             const = cst.ProcConst(nodes, MG, const)
             inipar[16,0] = 1
         else:     
             KG, MG, CG = cst.KG_MG_CG_condens(KG, MG, CG, const)
         # End if    
    # End if
    #
    # Compute internal forces considering the effective stiffness matrix, KE
    if DYN_flag == 1:
       KE  = effective(KG , MG , CG , ac)
       IGF = np.dot(KE, Up)
    else:
       IGF = np.dot(KG, Up)  
    # End if    
    #
    return KG, MG, CG, Msvar, IGF, ILF

def effective(KG , MG , CG , ac):
     
    KE = ac[3]*MG + ac[2]*CG + KG
    
    return KE     

def nodal_statloads(IBC, Nodal_loads, neq, ninc, inipar):
    """Assembles the global Right Hand Side Vector RHS due to nodal static loads
    
    Parameters
    ----------
    IBC : ndarray (int)
      Array that maps the nodes with number of equations.
    MG  : ndarray (floats)
      Global mass matrix of the system
    Nodal_loads : ndarray (floats)
      Static nodal loads
    neq : int
      Number of equations in the system after removing the nodes
      with imposed displacements.
    ninc : int
      Number of time steps for the solution
    inipar: ndarray (floats)  
      Analysis' Initial global parameters    
      
    Returns
    -------
    RHS : ndarray
      Array with the right hand side vector.
    
    """
    nloads = Nodal_loads.shape[0]
    RHS = np.zeros((neq, ninc))
    #
    NLSTA = int(inipar[5,0])
    #
    for i in range(nloads):
        node =  int(Nodal_loads[i,0])
        for j in range (6):
            IBC_load = IBC[node,j]
            if IBC_load != -1:
               if NLSTA == 0:  
                  RHS[IBC_load,:] = Nodal_loads[i,j+1]
               else:
                  load = np.linspace(0, Nodal_loads[i,j+1],ninc)  
                  RHS[IBC_load,:] = load[:]
               # End if
            # End if
        # End for j
    # End for i    
    #
    return RHS

def seismic_loadasem(IBC, MG, Seismo_signal, neq, ninc, inipar):
    """Assembles the global Right Hand Side Vector RHSG due to a seismic ground
       acceleration signal

    Parameters
    ----------
    IBC : ndarray (int)
      Array that maps the nodes with number of equations.
    MG  : ndarray (floats)
      Global mass matrix of the system
    Seismo_signal : ndarray (floats)
      Ground acceleration signal
    neq : int
      Number of equations in the system after removing the nodes
      with imposed displacements.
    ninc : int
      Number of time steps for the solution
    inipar: ndarray (floats)  
      Analysis' Initial global parameters  
      
    Returns
    -------
    RHSG : ndarray
      Array with the right hand side vector.
    """
    MG_exc   = np.array(MG)
    nnodes   = IBC.shape[0]
    ninc_sig = len(Seismo_signal)
    RHSG     = np.zeros((neq, ninc))
    #
    Accel_conditions = inipar[7:,0]
    aux_zeros        = np.zeros(neq)
    aux_ones         = np.ones(neq)
    #
    for i in range(6):
        for j in range (nnodes):
            #
            if IBC[j,i] != -1:
               #
               if Accel_conditions[i] == 0:
                  #  
                  MG_exc[IBC[j,i],:] = aux_zeros[:]
                  MG_exc[:,IBC[j,i]] = aux_zeros[:]
               # End if
            #End if 
        # End for j
    # End for i
    #
    m1 = np.dot(MG_exc,aux_ones)*-1
    flag = 0
    #
    for k in range(neq):
        #
        if ninc_sig <= ninc:
           RHSG[k, :ninc_sig] = m1[k]*Seismo_signal[:]
        else:
           RHSG[k, :] = m1[k]*Seismo_signal[:ninc]
           flag = 1
        # End if   
    # End for 
    if flag == 1:
         print(' ***** Warning: Total time of seismo signal is greater than solution total time *****')
    # End if
    #      
    return RHSG


def NewSignal(Tmin, DeltaT, Seismo_signal):
     
    """Computes a new signal data by linear interpolation
         
    Parameters
    ----------
    Tmin : float
      Minimum natural period of the system
    DeltaT: float
      Original time step
     Seismo_signal : ndarray (floats)
      Ground acceleration signal
      
    Returns
    -------
    New_Signal : ndarray
      Array with the updated signal data.
    """
    
    Ntimes_inter   = int(np.ceil(DeltaT/Tmin))      # Number of needed delta T in the actual interval
    New_DeltaT     = DeltaT/Ntimes_inter            # New time step
    Ntimes         = len(Seismo_signal)             # Number of time steps
    T_orig         = np.linspace(0, DeltaT*(Ntimes-1), Ntimes)      
    New_signal     = np.zeros(int((Ntimes-1)*Ntimes_inter + 1))
    New_signal[-1] = Seismo_signal[-1]
    #
    tcont = 0
    for i in range (Ntimes-1):
        m_i = (Seismo_signal[i+1] - Seismo_signal[i]) / (T_orig[i+1] - T_orig[i])
        Y1  = Seismo_signal[i]
        X1  = T_orig[i]
        #
        for j in range(Ntimes_inter):
            X = tcont*New_DeltaT
            Y = m_i*(X - X1) + Y1 
            New_signal[tcont] = Y
            tcont = tcont + 1
        # End for j
    # End for i  
    #
    return New_signal, New_DeltaT



def loadasem(IBC, MG, Seismo_signal, Nodal_loads, neq, ninc, inipar, Tmin, const):
    """Computes the global Right Hand Side Vector RHSG due to a seismic ground
       acceleration signal plus nodal static loads
    
    Parameters
    ----------
    IBC : ndarray (int)
      Array that maps the nodes with number of equations.
    MG  : ndarray (floats)
      Global mass matrix of the system
    Seismo_signal : ndarray (floats)
      Ground acceleration signal
    Nodal_loads : ndarray (floats)
      Static nodal loads  
    neq : int
      Number of equations in the system after removing the nodes
      with imposed displacements.
    ninc : int
      Number of time steps for the solution
    inipar: ndarray (floats)  
      Analysis' Initial global parameters     
    Tmin : int
      Minimum natural period of the system
    const: List
      Nodal constraints 
               
    Returns
    -------
    RHSG : ndarray
      Array with the right hand side vector.
    """
    
    NLSTA  = int(inipar[5,0])
    NLDYNA = int(inipar[6,0])
    DeltaT = inipar[0,0]
    DeltaT_fin = DeltaT
      

    if (NLSTA == 1) and (NLDYNA == 0):
         RHS_S  = nodal_statloads(IBC, Nodal_loads, neq, ninc, inipar)
         RHSG   = RHS_S
         #
    elif (NLSTA == 0) and (NLDYNA == 1):
         if Tmin < DeltaT:
              Seismo_signal, DeltaT_New = NewSignal(Tmin, DeltaT, Seismo_signal)
              print(' ***** Time step and seismo signal has been updated *****')
              DeltaT_fin  = DeltaT_New
              inipar[0,0] = DeltaT_fin
              ninc = int(inipar[1,0]/inipar[0,0])
         # End if
    # End if  
         RHS_S  = nodal_statloads(IBC, Nodal_loads, neq, ninc, inipar)
         RHS_D  = seismic_loadasem(IBC, MG, Seismo_signal, neq, ninc, inipar)
         RHSG   = RHS_S + RHS_D
         #
    elif (NLSTA == 1) and (NLDYNA == 1):
         inipar[5,0] = 0
         #
         if Tmin < DeltaT:
              Seismo_signal, DeltaT_New = NewSignal(Tmin, DeltaT, Seismo_signal)
              print(' ***** Time step and seismo signal has been updated *****')
              DeltaT_fin  = DeltaT_New
              inipar[0,0] = DeltaT_fin
              ninc = int(inipar[1,0]/inipar[0,0])
         # End if
         #
         RHS_S  = nodal_statloads(IBC, Nodal_loads, neq, ninc, inipar)
         RHS_D  = seismic_loadasem(IBC, MG, Seismo_signal, neq, ninc, inipar)
         RHSG   = RHS_S + RHS_D
         #
    # End if 
    #
    # Global condensation of RHSG
    #
    DPH_flag = const[0]
    CST_flag = const[1]
    #
    if (DPH_flag != 0) or (CST_flag != 0):
       RHSG = cst.RHSG_condens(RHSG, const)
    #End if   
    return RHSG, DeltaT_fin, ninc, inipar




