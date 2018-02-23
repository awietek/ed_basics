Exact Diagonalization Basics
=============================
.. image:: https://img.shields.io/badge/Python-2.7-brightgreen.svg?style=for-the-badge
.. image:: https://img.shields.io/badge/requires-NumPy%2C%20SciPy-blue.svg?style=for-the-badge

Tutorial for the FOR1807  Winter School 2018 in Marburg on
basic techniques in Exact Diagonalization. Basic codes are provided
that ought to be extended, as explained in the exercise sheet.

:author: Alexander Wietek
:license: GNU GPLv3
.. _hamiltonian_tfi.py: hamiltonian_tfi.py
.. _hamiltonian_hb_staggered.py: hamiltonian_hb_staggered.py
.. _hamiltonian_hb_xxz.py: hamiltonian_hb_xxz.py
.. _algorithm_lanczos.py: algorithm_lanczos.py
.. _example_groundstate_energy.py: example_groundstate_energy.py
.. _example_lanczos_convergence.py: example_lanczos_convergence.py
.. _exercises.pdf: exercises/exercises.pdf

Overview
-------------
This repository contains:

- Three files to create Hamiltonian matrices for three types of models:
  
  1. The Transverse Field Ising model on a linear chain,
     no Sz conservation or momentum conservation is used,
     hamiltonian_tfi.py_.
  2. The Heisenberg spin-1/2 model with staggered magnetic field
     on a linear chain, Sz conservation is used,
     hamiltonian_hb_staggered.py_.
  3. The Heisenberg XXZ model on a linear chain, Sz conservation and
     momentum conservation are used, hamiltonian_hb_xxz.py_.

- A basic implementation of the Lanczos algorithm, creating the T-matrix,
  algorithm_lanczos.py_.
  
- example scripts using the above functions example_groundstate_energy.py_
  and example_lanczos_convergence.py_.

- A set of exercises to be worked on during the hands-on session can be
  found in the file exercises.pdf_.

Installation
-------------
No installation is required. To use the python scripts, download
or clone into the git repository using

.. code-block:: bash
		
    $ git clone https://github.com/alexwie/ed_basics.git


Please *check in advance*, whether the scripts example_groundstate_energy.py_
and example_lanczos_convergence.py_ are working on your system:

.. code-block:: bash
		
    $ python example_groundstate_energy.py
    $ python example_lanczos_convergence.py
    
The numpy and scipy libraries are required. Please check, whether
you can update them to the latest version. If you are having any
problems please contact me.

Exercises
-------------
The exercises can be found in the file exercises.pdf_.

Developer
-------------
Alexander Wietek <alexander.wietek@uibk.ac.at>
