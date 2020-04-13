# -*- coding: utf-8 -*-
"""
Created by Juan Gomez and Julian Parra.

"""
from __future__ import division, print_function
from datetime import datetime
import preprocesor as pre
import assemutil as ass
import solutil as sol
import modalutil as modal
import constutil as cst
import postprocesor as pos

def Struct_DYN(folder):
    """  
    Parameters:
    ----------
    folder = Location of input files
    """
    #--------------------------------------------------------------------------------------------------------------------------------------
    # Pre-processing
    inipar, nodes, mats, elements, Nodal_loads, Msvar, ILF, Seismo_signal, const = pre.readin(folder)
    ninc, T, dt, ac, theta = pre.intparams(inipar)
    DME, IBC, neq = ass.DME(nodes, elements)
    const = cst.IBC_const(IBC, const)
    #--------------------------------------------------------------------------------------------------------------------------------------
    # System assembly
    Up, Vp                      = sol.inicond_U_V(neq)
    KG, MG, CG, Msvar, IGF, ILF = ass.assembler(inipar, Up, Msvar, ILF, elements, mats, nodes, neq, DME, ac, const)
    cst_neq, neq                = cst.neq_cst(const,neq)
    Up, Vp                      = sol.inicond_U_V(cst_neq) # New initial conditions after constraints
    Tfund, Tmin                 = modal.eigen(inipar, MG, KG)
    RHSG, DeltaT, ninc, inipar  = ass.loadasem(IBC, MG, Seismo_signal, Nodal_loads, neq, ninc, inipar, Tmin, const)
    ninc, T, dt, ac, theta      = pre.intparams(inipar)
    KG, MG, CG, Msvar, IGF, ILF = ass.assembler(inipar, Up, Msvar, ILF, elements, mats, nodes, neq, DME, ac, const)
    print("----------------------")
    print("Number of nodes: {}".format(nodes.shape[0]))
    print("Number of elements: {}".format(len(elements)))
    print("Number of equations: {}".format(neq))
    print("Number of equations after constraints: {}".format(cst_neq))
    print("----------------------")
    print("Natural periods of the system : ",Tfund)
    print("----------------------")
    print("Time step for solution: {} sec".format(DeltaT))
    print("Number of time increments: {}".format(ninc))
    print("----------------------")
    #--------------------------------------------------------------------------------------------------------------------------------------
    # System solution
    start_time = datetime.now()
    U, MvarsGen, ILFGen = sol.system_sol(inipar, Up, Vp, neq, RHSG, MG, KG, CG, DME, Msvar, ILF, elements, mats, nodes, ac, const, cst_neq)
    end_time   = datetime.now()
    print("Duration for the system's solution: {}".format(end_time - start_time))
    #--------------------------------------------------------------------------------------------------------------------------------------
    # Post-processing
    start_time = datetime.now()
    U = cst.Posproc_disp_cst(neq, ninc, const, U)
    end_time   = datetime.now()
    print('Duration for post processing: {}'.format(end_time - start_time))
    print("----------------------")
    print('Analysis terminated successfully!')
    print("----------------------")
    return (U, folder, IBC, nodes, elements, ninc, T, MvarsGen, ILFGen)
#------------------------------------------------------------------------------------------------------------------------------------------
# Execute
if __name__ == '__main__':
    displacement , folder , IBC , nodes, elements, ninc , T, MvarsGen, ILFGen = Struct_DYN('01_INPUT/')
    fig = pos.NodalDispPLT(displacement[0,:], T, ninc, ylabel = "Displacement")
    model = pos.GrafModel3D(elements, nodes)
    #td = pos.writeVTKs(elements,nodes,IBC,displacement)
    #histe = pos.PlasModel(MvarsGen, Element = 0, xlabel = "Y", ylabel = "P")  
