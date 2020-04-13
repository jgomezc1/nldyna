# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 06:36:40 2020

@author: JULIAN PARRA
"""
from __future__ import division, print_function
import numpy as np


def IBC_const(IBC, const):
    """
    This function gets the boundary conditions of the constrainted nodes and appends them to the constraints
    general information list.
    
    Parameters
    ----------
    IBC    : Boundary conditions with ID equation number
    const  : List with general information of the constraints in the assembly
    
    Output:
    -------
    const  : List with general information of the constraints in the assembly updated with the IBC.  
    
    """
    NCST    = len(const)-2
    IBC_CST = np.zeros((NCST,7))
    #
    for i in range (NCST):
        ID  = int(const[i+2][0])
        IBC_CST[i,0]  = const[i+2][0]
        IBC_CST[i,1:] = IBC[ID,:]
    # End for
    const.insert(2,IBC_CST)
    #      
    return const  

def ProcConst(nodes, MG, const):
    """
    This function:
         * Computes for diaphragms the center of mass X,Y,Z coordinates.
         * Relative coordinates to the master node of the restricted DOF.
         * Computes IBC for each master node
         * Computes ConstIBC that contains the IBC conditions for each node associated with
           each constraint
         
    Parameters:
    -----------
    nodes = array with nodal coordinates
    MG    = Global mass matrix
    const = list with type of constraints, IBC and DOF restrictions
     
    Output:
    -------
    ICC      = Index constraints conditions (0 = Diaphragms, 1 = General constraint)
    CMN      = array with master nodes coordinates (for diaphragms center of mass X,Y,Z coordinates).
    RLTC     = list of arrays with the relative coordinates to the master node of each constraint
    IBCc     = IBC conditions for each master node
    ConstIBC = IBC conditions for each node of each constraint
    
    """
    IBC_copy = np.array(const[2])
    CONS = np.array(const[3:])
    #
    if len(CONS) > 0:
       NDPH = int(max(CONS[:,1])+1)     #Total numer of diaphragms
       NCTS = int(abs(min(CONS[:,1])))  #Total number of constraints different of diaphragms
       CMN  = np.zeros((NDPH+NCTS,10))
       wc   = 0
       ICC  = []
       ConstIBC = [] # Storage of IBC conditions for each constraint
       IBC_copy = IBC_copy[np.argsort(CONS[:,1])]  # Sorting IBC by constraint ID
       CONS     = CONS[np.argsort(CONS[:,1])]      # Sorting CONS by constraint ID
    else:
       NDPH = 0
       NCTS = 0
       ICC  = []
       CMN  = []
       ConstIBC = []
       IBC_copy = []
       CONS = []
    # End if   
    RLTC = [] # List with relative distante to the master node
    #
    wc2 = 0   # Only for j cycle
    #
    for i in range (NDPH+NCTS):
        # 
        CMN[i,0] = CMN[i,0] + int(max(nodes[:,0])) + 1 + i
        IDcons   = CONS[i+wc,1]
        cn       = 0     # Number of  one constraint associatted nodes
        RLTCi    = []    # List with relative distance to the master node
        IBCci    = []    # List with IBC conditions for each constraint
        #
        if IDcons >= 0: #If the constraint is a diaphragm
           #
           CXM = 0     # X coord. * mass accumulator
           CYM = 0     # Y coord. * mass accumulator
           CZM = 0     # Z coord. * mass accumulator
           M   = 0     # Mass accumulator
           ICC.append(0)
           #
           while (CONS[wc,1] == IDcons):
                 cn = cn + 1
                 IBCci.append(IBC_copy[wc])
                 IDn  = int(CONS[wc,0])
                 IBCi = int(max(IBC_copy[wc])) #Should be max because one can have restricted DOF (IBC = -1)
                 #
                 CXM = CXM + nodes[IDn,1]*MG[IBCi,IBCi]
                 CYM = CYM + nodes[IDn,2]*MG[IBCi,IBCi]
                 CZM = CZM + nodes[IDn,3]*MG[IBCi,IBCi]
                 M   = M   + MG[IBCi,IBCi]
                 #
                 if wc == (len(CONS)-1):
                    wc = wc + 1
                    break
                 # End if
                 wc = wc + 1
           # End while
           #
           CMN[i,1:4] = [CXM/M, CYM/M, CZM/M] # Coordinates
           CMN[i,4:]  = CONS[wc-1,3:]
        else:
           ICC.append(1)  
           while (CONS[wc,1] == IDcons):
                 cn = cn + 1
                 IBCci.append(IBC_copy[wc])
                 IDn  = int(CONS[wc,0])
                 if CONS[wc,2] == 1: # If it is the master node
                    CMN[i,1:4] = nodes[IDn,1:4]
                    CMN[i,4:]  = CONS[wc,3:]
                 # End if
                 #
                 if wc == (len(CONS)-1):
                    wc = wc + 1
                    break
                 # End if
                 wc = wc + 1
           # End while
        # End if
        #
        ConstIBC.append(np.array(IBCci))
        # Compute relative coordinates to the master node
        for j in range (cn):
            #
            IDn = int(CONS[j+wc2,0])
            RCCMi = [] #Relative coordinates to the master node 
            RCCMi.append(nodes[IDn,0]) 
            RCCMi.append(nodes[IDn,1] - CMN[i,1]) # X coord. of node to the X cord of master node
            RCCMi.append(nodes[IDn,2] - CMN[i,2]) # Y coord. of node to the Y cord of master node
            RCCMi.append(nodes[IDn,3] - CMN[i,3]) # Z coord. of node to the Z cord of master node
            RLTCi.append(RCCMi)
        # End for j
        wc2   = wc
        RLTCi = np.array(RLTCi)
        RLTC.append(RLTCi)
    #End for i
    #
    # Compute index boundary conditions for each master node
    NDOF   = len(MG)
    IBCc   = np.zeros((len(CMN),7))
    IBCc   = np.where(IBCc==0,-1,IBCc)
    ic     = 0
    for i in range (len(CMN)):
        IBCc[i,0] = CMN[i,0] 
        #
        for j in range (6):
            if CMN[i,j+4] != 0:
               IBCc[i,j+1] = NDOF + ic
               ic = ic + 1
            # End if
        # End for j
    # End for i
    #
    del(const[3:])
    const.append(CONS)
    const.append(ICC)
    const.append(CMN)
    const.append(RLTC) 
    const.append(IBCc) 
    const.append(ConstIBC)
    #
    # ID of constraints DOF for condensation, this will be used for deleting rows and columns of
    # KG, MG, CG, RHSG matrixes.
    #
    DOF_delete = []
    for i in range (len(IBCc)):
        IBCci      = IBCc[i]
        NodesIBCci = ConstIBC[i]
        # 
        for j in range (6):
            # 
            if IBCci[j+1] != -1:
               #
               for k in range (len(NodesIBCci)):
                   IDdof = int(NodesIBCci[k,j+1])
                   DOF_delete.append(IDdof)
               # End for k
             # End if   
        # End for j
    # End for i
    DOF_stay = np.linspace(0,len(MG)-1,len(MG),dtype=int)
    DOF_stay = np.delete(DOF_stay, DOF_delete, 0)
    DOF_proc = []
    DOF_proc.append(np.array(DOF_delete))
    DOF_proc.append(DOF_stay)
    const.append(DOF_proc)        
    #   
    return const           

def KG_MG_CG_condens(KG, MG, CG, const):
    """
    This function computes condensed KG, MG, DG matrixes.
    
    Parameters:
    -----------
    KG    = Global stiffness matrix
    MG    = Global mass matrix
    CG    = Global damping matrix
    const = Global list with general processing of constraint information:
         CMN   = array with master nodes coordinates and boundary DOF restrictions
         IBCc  = IBC conditions for each master node
         ConstIBC = IBC conditions for each node of each constraint
     
    Output:
    -------
    KGc, MGc, CGc = Condensed stiffness, mass and damping matrixes.
     
    """
    #
    if (const[0] != 0) or (const[1] != 0): #If there is any constraint
        CMN = const[5]
        IBCc = const[7]
        ConstIBC = const[8]
        DOF_delete = const[9][0]
        #
        NDOFc  = int(np.sum(CMN[:,4:]))
        NDOF   = len(KG)
        KGc    = np.zeros((NDOF+NDOFc,NDOF+NDOFc))
        MGc    = np.zeros_like(KGc)
        CGc    = np.zeros_like(KGc)
        #
        KGc[:NDOF,:NDOF] = KG[:,:]
        MGc[:NDOF,:NDOF] = MG[:,:]
        CGc[:NDOF,:NDOF] = CG[:,:]
        #
        for i in range (len(IBCc)):
            IBCci      = IBCc[i]
            NodesIBCci = ConstIBC[i]
            # 
            for j in range (6):
                # 
                if IBCci[j+1] != -1:
                   IDdof_c = int(IBCci[j+1])
                   #
                   for k in range (len(NodesIBCci)):
                       IDdof = int(NodesIBCci[k,j+1])
                       #
                       KGc[IDdof_c,:] = KGc[IDdof_c,:] + KGc[IDdof,:]
                       KGc[:,IDdof_c] = KGc[:,IDdof_c] + KGc[:,IDdof]
                       #
                       MGc[IDdof_c,:] = MGc[IDdof_c,:] + MGc[IDdof,:]
                       MGc[:,IDdof_c] = MGc[:,IDdof_c] + MGc[:,IDdof]
                       #
                       CGc[IDdof_c,:] = CGc[IDdof_c,:] + CGc[IDdof,:]
                       CGc[:,IDdof_c] = CGc[:,IDdof_c] + CGc[:,IDdof]
                       #
                  # End for k
                # End if
            #End for j
        # End for i
        #
        KGc = np.delete(KGc, DOF_delete, 0)
        KGc = np.delete(KGc, DOF_delete, 1)
        #
        MGc = np.delete(MGc, DOF_delete, 0)
        MGc = np.delete(MGc, DOF_delete, 1)
        #
        CGc = np.delete(CGc, DOF_delete, 0)
        CGc = np.delete(CGc, DOF_delete, 1) 
        #
    else:
        KGc = KG
        MGc = MG
        CGc = CG
    # End if
    #              
    return KGc, MGc, CGc     

def RHSG_condens(RHSG, const):
    """
    This function computes condensed RHSG (loads vector for each deltaT).
    
    Parameters:
    -----------
    RHSG  = Global nodal loads
    const = Global list with general processing of constraint information:
         CMN   = array with master nodes coordinates and boundary DOF restrictions
         IBCc  = IBC conditions for each master node
         ConstIBC = IBC conditions for each node of each constraint
     
    Output:
    -------
    RHSGc = Condensed RHSG vector (loads vector for each deltaT).
     
    """
    #
    if (const[0] != 0) or (const[1] != 0): #If there is any constraint
       CMN      = const[5]
       RLTC     = const[6]
       IBCc     = const[7]
       ConstIBC = const[8]
       DOF_delete = const[9][0]
       #
       NDOFc  = int(np.sum(CMN[:,4:]))
       NDOF   = len(RHSG)
       RHSGc  = np.zeros((NDOF+NDOFc,len(RHSG[0])))
       #
       RHSGc[:NDOF,:] = RHSG[:,:]
       #
       for i in range (len(IBCc)):
           IBCci      = IBCc[i]
           NodesIBCci = ConstIBC[i]
           RLTCci     = RLTC[i]
           #
           Fx_dy = np.zeros(len(RHSG[0])) # -Mzz
           Fx_dz = np.zeros(len(RHSG[0])) #  Myy
           #
           Fy_dx = np.zeros(len(RHSG[0])) #  Mzz
           Fy_dz = np.zeros(len(RHSG[0])) # -Mxx 
           #
           Fz_dy = np.zeros(len(RHSG[0])) #  Mxx
           Fz_dx = np.zeros(len(RHSG[0])) # -Myy              
           #
           for j in range (6):
               # 
               if IBCci[j+1] != -1:
                  IDdof_c = int(IBCci[j+1])
                  #
                  for k in range (len(NodesIBCci)):
                      #
                      IDdof = int(NodesIBCci[k,j+1])
                      # 
                      if (j <= 2): # Forces
                          # 
                          RHSGc[IDdof_c,:] = RHSGc[IDdof_c,:] + RHSGc[IDdof,:]
                          if j == 0:   # Fx
                             Fx_dy[:] = Fx_dy[:] - RHSGc[IDdof,:]*RLTCci[k,2]
                             Fx_dz[:] = Fx_dz[:] + RHSGc[IDdof,:]*RLTCci[k,3]
                          elif j == 1: # Fy
                             Fy_dx[:] = Fy_dx[:] + RHSGc[IDdof,:]*RLTCci[k,1]
                             Fy_dz[:] = Fy_dz[:] - RHSGc[IDdof,:]*RLTCci[k,3]
                          else:
                             Fz_dy[:] = Fz_dy[:] + RHSGc[IDdof,:]*RLTCci[k,2]
                             Fz_dx[:] = Fz_dx[:] - RHSGc[IDdof,:]*RLTCci[k,1]
                          # End if   
                          #
                      else:        # Moments 
                          #
                          RHSGc[IDdof_c,:] = RHSGc[IDdof_c,:] + RHSGc[IDdof,:] 
                      # End if
                      #
                  # End for k
                  if   j == 3:
                       RHSGc[IDdof_c,:] = RHSGc[IDdof_c,:] + Fz_dy + Fy_dz
                  elif j == 4:
                       RHSGc[IDdof_c,:] = RHSGc[IDdof_c,:] + Fx_dz + Fz_dx
                  elif j == 5:
                       RHSGc[IDdof_c,:] = RHSGc[IDdof_c,:] + Fy_dx + Fx_dy
                  # End if
                  #
               # End if
           #End for j
       # End for i
       RHSGc = np.delete(RHSGc, DOF_delete, 0)
       #
    else:
       RHSGc = RHSG
    # End if 
    #
    return RHSGc


def disp_retrieve(neq, const, Up):
    """
    This function computes the element's displacements after matrixes condensation,
    it retrieves nodal displacements after the contraint system solution.
    
    Parameters:
    -----------
    neq   = Number of system's equations (without constraints).
    const = Global list with general processing of constraint information.
    Up    = Global displacements considering constraints.
     
    Output:
    -------
    Upn = New displacements for each element (after constraints).
    
    """

    Upn   = np.zeros((neq,2))
    IBCcn = np.array(const[7])

    DeltaDOF = max(IBCcn[-1,1:]) - (len(Up)-1)

    IBCcn[:,1:] = IBCcn[:,1:] - DeltaDOF
    IBCcn       = np.where(IBCcn<0,-1,IBCcn)
    NCST        = len(const[7])
    DOF_stay    = const[9][1]
    #
    Upn[:len(DOF_stay),1] = Up[:len(DOF_stay)]
    Upn[:len(DOF_stay),0] = DOF_stay[:]
    ConstIBC              = const[8]
    RLTC                  = const[6]
    #
    cont                  = 0
    #
    for i in range (NCST):
        IBCcn_i   = IBCcn[i,:]
        ConstIBCi = ConstIBC[i]
        RLTCci    = RLTC[i]
        #
        for j in range (len(ConstIBCi)):
            dx = RLTCci[j][1]
            dy = RLTCci[j][2]
            dz = RLTCci[j][3]
            #
            for k in range (6):
                #
                if IBCcn_i[k+1] != -1:
                   #  
                   IDdof_k   = int(ConstIBCi[j,k+1])
                   IDdofcn_i = int(IBCcn_i[k+1])
                   #
                   Upn[len(DOF_stay)+cont,0] = IDdof_k
                   #
                   if k == 0:  # Traslation along X axis
                      tx = Up[IDdofcn_i]
                      #
                      if IBCcn_i[6] != -1:
                         rz = Up[int(IBCcn_i[6])]
                      else:
                         rz = 0
                      # End if
                      if IBCcn_i[5] != -1:
                         ry = Up[int(IBCcn_i[5])]
                      else:
                         ry = 0
                      # End if
                      Upn[len(DOF_stay)+cont,1] = tx - rz*dy + ry*dz
                      #
                   elif k == 1:  # Traslation along Y axis
                      ty = Up[IDdofcn_i]
                      #
                      if IBCcn_i[6] != -1:
                         rz = Up[int(IBCcn_i[6])]
                      else:
                         rz = 0
                      # End if
                      if IBCcn_i[4] != -1:
                         rx = Up[int(IBCcn_i[4])]
                      else:
                         rx = 0
                      # End if
                      Upn[len(DOF_stay)+cont,1] = ty + rz*dx - rx*dz
                      #
                   elif k == 2: # Traslation along Z axis
                      tz = Up[IDdofcn_i]
                      #
                      if IBCcn_i[4] != -1:
                         rx = Up[int(IBCcn_i[4])]
                      else:
                         rx = 0
                      # End if
                      if IBCcn_i[5] != -1:
                         ry = Up[int(IBCcn_i[5])]
                      else:
                         ry = 0
                      # End if
                      Upn[len(DOF_stay)+cont,1] = tz + rx*dy - ry*dz
                      #
                   else:
                      Upn[len(DOF_stay)+cont,1] = Up[IDdofcn_i]
                   #End if              
                   #
                   cont = cont + 1
                   #
    Upn = Upn[np.argsort(Upn[:,0])]  # Sorting IBC by constraint ID
    Upn = np.delete(Upn, 0, axis=1)
    #
    return Upn 

def neq_cst(const,neq):
    """
    This function computes the system's new number of equations after constraints.
    
    Parameters:
    -----------
    neq   = Number of system's equations (without constraints).
    const = Global list with general processing of constraint information.
     
    Output:
    -------
    neq = New number of equations
    
    """
    #
    DPH_flag = const[0]
    CST_flag = const[1]
    #
    if (DPH_flag != 0) or (CST_flag != 0):
       initial_IBC = const[5]
       neq_stay = len(const[9][1])
       neq_const = np.sum(initial_IBC[:,4:])
       new_neq = int(neq_stay + neq_const)
    else:
       new_neq = neq
    # End if   
    #
    return  new_neq, neq

def Posproc_disp_cst(neq, ninc, const, U):
    """
    This function computes general nodal displacements for each time increment
    
    Parameters:
    -----------
    neq   = Number of system's equations (without constraints).
    ninc  = Number of time steps in system's solution
    const = Global list with general processing of constraint information.
    U     = list with global DOF displacements
     
    Output:
    -------
    Un = New array with general nodal displacements
    
    """
    DPH_flag = const[0]
    CST_flag = const[1]
    #
    if (DPH_flag != 0) or (CST_flag != 0):
         #
       Un= np.zeros((neq,ninc))
       #
       for i in range (ninc):
           Upi = U[:,i]
           Upi = disp_retrieve(neq, const, Upi)
           Un[:,i] = Upi[:,0]
       # End for
    else:
       Un = U  
    # End if
    #
    return Un    
    

    
     

      
          


            







    
    
    
    