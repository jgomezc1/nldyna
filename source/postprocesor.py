# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 06:22:56 2018

@author: JULIAN PARRA
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Post-processing functions for graphics, and result visualization

def NodalDispPLT(DOF,TotalT,ninc,ylabel):
    
    """ This function plots displacement Vs. t
    INPUT:
    -----
    - DOF   : Displacement response history (1d array) of the degree of freedom
    - TotalT: Total time of the step by step integration procedure of solution
    
    OUTPUT:
    ------
    - python plot displacement Vs. t
    
    """

    t = np.linspace(0,TotalT,len(DOF))
    plt.figure(figsize=(6.7,4))
    plt.plot(t,DOF,'gray')
    plt.grid(True)
    plt.xlim(xmin=0,xmax=TotalT)
    plt.title("Displacement history for the specified DOF")
    plt.xlabel("Time (sec)")
    plt.ylabel(ylabel)
    plt.show()
    return
    
def GrafModel(Elements,Nodes):
    
    """ This function plots Model, only for frame elements
    INPUT:
    -----
    - Elements: Element conectivity (array)
    - Nodes: Nodes coordinates (array)

    OUTPUT:
    ------
    - python model plot
    
    """
    Nlines = len(Elements)
    #
    plt.figure(figsize=(7,4))
    for i in range (Nlines):
        Cordx = np.array([Nodes[Elements[i][3]][1],Nodes[Elements[i][4]][1]])
        Cordy = np.array([Nodes[Elements[i][3]][2],Nodes[Elements[i][4]][2]])
        plt.plot(Cordx,Cordy,'black')
    #End for
    plt.xlim(min(Nodes[:,1])-1,max(Nodes[:,1])+1)
    plt.ylim(min(Nodes[:,2])-1,max(Nodes[:,2])+1)
    plt.xlabel("Y")
    plt.ylabel("X")
    # 
    plt.show()
    return

def GrafModel3D(Elements,Nodes):
    
    """ This function plots 3D Model, only for frame elements
    INPUT:
    -----
    - Elements: Element conectivity (array)
    - Nodes: Nodes coordinates (array)

    OUTPUT:
    ------
    - python model plot
    
    """
    Nlines = len(Elements)
    #
    fig = plt.figure(figsize=(7,4))
    ax  = fig.gca(projection = '3d')
    for i in range (Nlines):
        Cordx = np.array([Nodes[Elements[i][3]][1],Nodes[Elements[i][4]][1]])
        Cordy = np.array([Nodes[Elements[i][3]][2],Nodes[Elements[i][4]][2]])
        Cordz = np.array([Nodes[Elements[i][3]][3],Nodes[Elements[i][4]][3]])
        ax.plot(Cordx,Cordy,Cordz,'k')
    #End for
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    #
    ax.set_xlim3d(min(Nodes[:,1])-1,max(Nodes[:,1])+1)
    ax.set_ylim3d(min(Nodes[:,2])-1,max(Nodes[:,2])+1)
    ax.set_zlim3d(min(Nodes[:,3]),max(Nodes[:,3])+1)
    return

def PlasModel(MvarsGen, Element, xlabel, ylabel):
     
    """ This function plots from results the elasto-plastic histeretic curve
    INPUT:
    -----
    - MsvarGen: Python list. It storages the history of state variables of each element
    - Element : Integer. 
    - xlabel  : String for title of X axis
    - ylabel  : String for title of Y axis
    
    OUTPUT:
    ------
    - elastoplatic curve plot
    
    """
    X = np.zeros(len(MvarsGen))
    Y = np.zeros(len(MvarsGen))
    
    for i in range (len(MvarsGen)):
        Mvars = MvarsGen[i] 
        X[i] = Mvars[Element][1]
        Y[i] = Mvars[Element][0]
         
    plt.figure(figsize=(6.7,4))    
    plt.plot(X,Y,'gray')
    Title = "Constitutive model history for element "
    plt.title(Title + str(Element))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.show()

    return
    
def writeVTKs(Elements,Nodes,IBC,disp):
    
    """ This function save vtks files for model animation with paraview
    INPUT:
    -----
    - Elements: Element conectivity (array)
    - Nodes: Nodes coordinates (array)
    - IBC: Index boundary condition (array)
    - disp: displacement history of non restraint nodes (array)
    
    OUTPUT:
    ------
    - VTK files
    
    """
    #
    Ninc = len(disp[0,:])
    Nnodes = len(Nodes)
    Nelem = len(Elements)
    #
    Totaldisp =np.zeros((Nnodes*2,Ninc)) # For all degrees of freedom, even those with restraints
    #
    for i in range (Nnodes):
        for j in range (2):
            #
            if IBC[i][j] != -1:
                Totaldisp[2*i + j,:] = disp[IBC[i][j],:]
            #End if
        # End for j
    #End for i
    #
    # For each deltaT, its generated a VTK file with coordinates for each node acummulatting nodal displacements x and y
    #
    for i in range(Ninc):
        VTKi = open('03_VTKs/' + 't' + str(i) + '.vtk','w')
        VTKi.write('# vtk DataFile Version 2.0\n')
        VTKi.write('File for t = ' + str(i) + '\n')
        VTKi.write('ASCII\n')
        VTKi.write('DATASET UNSTRUCTURED_GRID\n')
        VTKi.write('POINTS ' + str(Nnodes) + ' float\n')
        #
        for k in range (Nnodes):
            VTKi.write('%10.2f %10.2f %10.2f\n' %(Nodes[k][1],Nodes[k][2],0.0))
        # End for k
        VTKi.write('\n')
        VTKi.write('CELLS ' + str(Nelem) + ' ' + str(Nelem*3) + '\n')
        #
        for k in range (Nelem):
            VTKi.write('%10i %10i %10i\n' %(2,Elements[k][3],Elements[k][4]))
        # End for k
        VTKi.write('\n')
        VTKi.write('CELL_TYPES ' + str(Nelem) + '\n')
        #
        for k in range (Nelem):
            VTKi.write('%10i\n' %(3))
        # End for k
        VTKi.write('\n')
        VTKi.write('POINT_DATA ' + str(Nnodes)+ '\n')
        VTKi.write('SCALARS dispX float\n')
        VTKi.write('LOOKUP_TABLE default\n')
        #
        FMT=1*'%10.2f'
        for k in range (Nnodes):
            VTKi.write(FMT %(Totaldisp[2*k][i]) + '\n')
        # End for k
        #
        VTKi.write('\n')
        VTKi.write('SCALARS dispY float\n')
        VTKi.write('LOOKUP_TABLE default\n')
        #
        FMT=1*'%10.2f'
        for k in range (Nnodes):
            VTKi.write(FMT %(Totaldisp[2*k + 1][i]) + '\n')
        # End for k
        VTKi.close
    # End for i

    return Totaldisp
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    