Exact Diagonalization Basics
=============================

Basic scripts for diagonalizing simple lattice models and extracting
physical information.

:author: Alexander Wietek
:license: GNU GPLv3

Overview
-------------
This repository contains

- A python package to create Hamiltonian matrices for three types of models
  1. The Transverse Field Ising model on a linear chain,
     no Sz conservation or momentum conservation is used
  2. The Heisenberg spin-1/2 model with staggered magnetic field
     on a linear chain, Sz conservation is used
  3. The Heisenberg XXZ model on a linear chain, Sz conservation and
     momentum conservation are used
  the package can be found in the subfolder hamiltonians
- A basic implementation of the Lanczos algorithm, creating the T-matrix,
  can be found in the folder algorithms
- example scripts using the above packages to evaluate properties of the
  models, these are contained in the folder examples and also shown in
  the section `examples`_
  .. _examples: examples/examples.rst
- exercises, to deepen the understanding of creating Hamiltonians and
  employing algorithms to evaluate physics

Installation
-------------
No installation is required. To use the python scripts, download
or clone into the git repository using

.. code-block:: bash
		
    $ git clone git@github.com:alexwie/ed_basics.git

Exercises
-------------
A collection of exercises can be found in the `exercise`_ section.

.. _exercise: exercises/exercises.rst

Developer
-------------
Alexander Wietek <alexander.wietek@uibk.ac.at>
