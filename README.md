# NLDYNA - Nonlinear Dynamic Analysis.
    
## By Julian David Parra Cardona.

## This document describes the results of the project conducted by the author in partial fulfillment of the requirements to complete the degree of Master of Engineering at Universidad EAFIT.

## Introduction

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jgomezc1/nldyna/master)

**(Click on the Binder badge to execute the NBs in an external Jupyter server)**

The study of structural systems considering inelastic response was, until very recently, only possible through the use of commercial packages with subsequent limitations as research tools. Commercial packages are very rigid since not all of them allow the addition of independent elements and/or material models as often required in research activities. At the same time the used solvers in commercial codes are of the black-box type making its extension to general tasks an impossible goal. The recent development of high level lenguages (e.g., Python) facilitates the development of highly efficient in-house implementations. This project describes a general in-house finite element assembler and solver aimed at studying the non-linear response of dynamic systems. The code is intended to be used in the testing of material models, complex kinematic formulations and novel structural systems commonly required in research activities. The code has the following features:

* **It is multiphysics oriented:** The code is a general dynamic Newton-Raphson solver where the physical context is provided by the user in terms of material and/or element models.

* **The implementation has been fully parametrized:** It does not have an implicit space dimensionality and problems with an arbitrary number of degrees of freedom per node can be solved.

* **Python based user elements and material models:** The implementation of user elements and user constitutive models is highly simplified in comparisson with commercial codes as it is conducted in a high level language.

* **Easily coupled with independent scripts:** Since the code is fully open and written in a modular structure it can be coupled with external independent scripts required in specific analysis and design problems.


## Nonlinear dynamic analysis of generalized finite element problems.
**NLDYNA** is a generalized finite element program for the solution of time-dependent non-linear problems. The code is able to handle static and dynamic analysis problems assumed of hyperbolic nature. It is generalized as it can solve user defined problems in different physical contexts through the implementation of user elements and user constitutive responses.

In **NLDYNA** a dynamic problem is splitted into several time increments and each increment is solved by a Newton-Raphson algorithm. Time stepping is conducted by an implicit Wilson $\theta$-method. The solution of linear static problems takes place in a single increment and a single iteration. 

A model is defined in **NLDYNA** through 5 easy to write input data files containing: (i) basic problem parameters (ii) nodal point data (iii) element data (iv) loads data and (v) material data. The model can use elements available in the code's own library, specific user defined elements or a combination of both. Similarly, a model can also use **NLDYNA's** available elements in combination with user defined material models.

![Shaking in 3D building.](./notebooks/img/Model_Page.png)

## Contents:

* [Formulation: Description of the Newton-Raphson incremetal scheme.](https://nbviewer.jupyter.org/github/jgomezc1/nldyna/blob/master/notebooks/02_Formulation.ipynb)

* [Input files: Template of input files required to create a model.](https://nbviewer.jupyter.org/github/jgomezc1/nldyna/blob/master/notebooks/03_NLDYNA.ipynb)

* [User elements: Example of a user defined element subroutine.](https://nbviewer.jupyter.org/github/jgomezc1/nldyna/blob/master/notebooks/04_UEL_subroutine.ipynb)

* [User materials: Example of a user defined material subroutine.](https://nbviewer.jupyter.org/github/jgomezc1/nldyna/blob/master/notebooks/05_UMAT_subroutine.ipynb)

## Examples

* [Static analysis of a 2D-frame: Basic linear analysis of a 2D frame.](https://nbviewer.jupyter.org/github/jgomezc1/nldyna/blob/master/notebooks/06_Example01.ipynb)

* [Dynamic analysis of a 2D-frame: Basic linear dynamic analysis of a 2D frame.](https://nbviewer.jupyter.org/github/jgomezc1/nldyna/blob/master/notebooks/07_Example02.ipynb)

* [Rigid diaphragm in 3D-frame: Dynamic analysis of a 3D frame with rigid diaprhagm.](https://nbviewer.jupyter.org/github/jgomezc1/nldyna/blob/master/notebooks/08_Example03.ipynb)

* [Non-linear spring: One-dimensional non-linear spring element.](https://nbviewer.jupyter.org/github/jgomezc1/nldyna/blob/master/notebooks/09_Example04.ipynb)

* [Based-isolated 2D frame: Dynamic analysis of based-isolated 2D frame](https://nbviewer.jupyter.org/github/jgomezc1/nldyna/blob/master/notebooks/10_Example05.ipynb)

* [Laterally loaded pile: Non-linear analysis of a laterally loaded pile foundation](https://nbviewer.jupyter.org/github/jgomezc1/nldyna/blob/master/notebooks/11_Example06.ipynb)

## SolidsPy.
This program was written on top of SolidsPy, the 2D Finite Element Analysis By Juan Gomez and Nicolas Guarín-Zapata.
https://github.com/AppliedMechanics-EAFIT/SolidsPy
