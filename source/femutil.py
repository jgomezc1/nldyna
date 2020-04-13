 # -*- coding: utf-8 -*-
"""
femutil.py
----------

"""
from __future__ import division, print_function

def eletype(iet):
    """Assigns number to degrees of freedom

    According to iet assigns number of degrees of freedom, number of
    nodes and minimum required number of integration points.

    Parameters
    ----------
    iet :  int
      Type of element. These are:
        0.  Linear - 1D Spring in X direction.   
        1.  Linear - 2D Simple frame.
        2.  Linear - 2D Full frame.
        3.  Linear - 2D Truss.
        4.  NonLin - 4 noded plate.
        5.  NonLin - 1D Spring in X direction.
        6.  Linear - 2D Shear-Rotational spring.
        7.  Linear - 1D Rotational spring.
        8.  NonLin - 1D Rotational spring.
        9.  NonLin - 1D pile spring P-Y curve for stiff soil
        10. Linear - 3D full frame bending  ans shear effects

    Returns
    -------
    ndof : int
      Number of degrees of freedom for the selected element.
    nnodes : int
      Number of nodes for the selected element.
    ngpts : int
      Number of Gauss points for the selected element.

    """
    if iet == 0:
        ndof = 2
        nnodes = 2
        ngpts = 1
    elif iet == 1:
        ndof = 6
        nnodes = 2
        ngpts = 1
    elif iet == 2:
        ndof = 6
        nnodes = 2
        ngpts = 1
    elif iet == 3:
        ndof = 4
        nnodes = 2
        ngpts = 1
    elif iet == 4:
        ndof = 8
        nnodes = 4
        ngpts = 4
    elif iet == 5:
        ndof = 2
        nnodes = 2
        ngpts = 1
    elif iet == 6:
        ndof = 4
        nnodes = 2
        ngpts = 1
    elif iet == 7:
        ndof = 2
        nnodes = 2
        ngpts = 1
    elif iet == 8:
        ndof = 2
        nnodes = 2
        ngpts = 1
    elif iet == 9:
        ndof = 2
        nnodes = 2
        ngpts = 1    
    elif iet == 10:
        ndof = 12
        nnodes = 2
        ngpts = 1  
    return ndof, nnodes, ngpts

if __name__ == "__main__":
    import doctest
    doctest.testmod()
