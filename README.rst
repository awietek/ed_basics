Exact Diagonalization Basics
=============================

Tutorial for the FOR1807  Winter School 2018 in Marburg on
basic techniques in Exact Diagonalization. Basic codes are provided
that ought to be extended, as explained in the exercise sheet.

:author: Alexander Wietek
:license: GNU GPLv3

Overview
-------------
This repository contains

- Three files to create Hamiltonian matrices for three types of models
  
  1. The Transverse Field Ising model on a linear chain,
     no Sz conservation or momentum conservation is used,
     `<hamiltonian_tfi.py>`
  2. The Heisenberg spin-1/2 model with staggered magnetic field
     on a linear chain, Sz conservation is used,
     `<hamiltonian_hb_staggered.py>`
  3. The Heisenberg XXZ model on a linear chain, Sz conservation and
     momentum conservation are used,
     `<hamiltonian_hb_xxz.py>`
     
- A basic implementation of the Lanczos algorithm, creating the T-matrix,
  `<algorithm_lanczos.py>`
  
- example scripts using the above functions `<example_groundstate_energy.py>`
  and `<example_lanczos_convergence.py>`
  
- A set of exercises to be worked on during the hands-on session, `<exercises/exercises.pdf>`

Installation
-------------
No installation is required. To use the python scripts, download
or clone into the git repository using

.. code-block:: bash
		
    $ git clone git@github.com:alexwie/ed_basics.git

Please check in advance, whether the scripts `<example_groundstate_energy.py>`
  and `<example_lanczos_convergence.py>` are working on your system.

Exercises
-------------
The exercises can be found on the manuscript 

Developer
-------------
Alexander Wietek <alexander.wietek@uibk.ac.at>
