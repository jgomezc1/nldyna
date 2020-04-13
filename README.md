# NLDYNA - Nonlinear Dynamic Analysis program.
## Introduction

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jgomezc1/nldyna/master)

The study of structural systems considering inelastic response was, until very recently, only possible through the use of commercial packages with subsequent limitations as a research tool. For instance, commercial packages are very rigid since not all of them allow the addition of independent elements and/or material models as often required in research activities. In top of that restriction the used solvers in commercial codes are of the black-box type making its extension to general tasks an impossible goal. With the recent development of high level lenguages (like Python) it is now possible to develop very efficient in-house implementations. This project describes a general in-house finite element assembler and solver aimed at studying the non-linear response of dynamic systems. The code is intended to be used in the testing of material models and/or complex kinematic formulations commonly required in research activities. The code has the following advantages:

* **It is multiphysics oriented:** The code is just a general dynamic Newton-Raphson solver where the physical context is provided by the user in terms of material and/or element models.

* **The implementation has been fully parametrized:** It does not have an implicit space dimensionality and problems with an arbitrary number of degrees of freedom per node can be solved.

* **Python based user elements and material models:** The implementation of user elements and user constitutive models is highly simplified in comparisson with commercial codes as it is connducted in a high level language like Python.

* **Easily coupled with independent scripts:** Since the code is fully open and written in a modular structure it can be coupled with external independent scripts required in specific analysis and design problems.


## Nonlinear dynamic analysis of generalized finite element problems
**NLDYNA** is a generalized finite element program for the solution of time-dependent non-linear problems. The code is able to handle static and dynamic analysis problems assumed of hyperbolic nature. It is generalized as it can solve user defined problems in different physical contexts through the implementation of user elements and user constitutive responses.

In **NLDYNA** a dynamic problem is splitted into several time increments and each increment is solved by a Newton-Raphson algorithm. Time stepping is conducted by an implicit Wilson $\theta$-method. The solution of linear static problems takes place in a single increment and a single iteration. 

A model is defined in **NLDYNA** through 6 easy to write input data files containing basic problem parameters, nodal point data, element data and loads data. The model can use elements available in the code's library, specific user defined elements or a combination of both. Similarly, a model can also use **NLDYNA's** available elements in combination with user defined material models.

![Shaking in 3D building.](./notebooks/img/Model_Page.png)
