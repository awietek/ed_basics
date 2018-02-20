# -*- coding: utf-8 -*-
'''
Example how to use the Lanczos algorithm in combination with the
sparse matrix data from creating Hamiltonians

:author: Alexander Wietek
:email: alexander.wietek@uibk.ac.at
:year: 2018
'''   

from __future__ import absolute_import, division, print_function

import numpy as np
import scipy as sp
from scipy import sparse

import hamiltonian_hb_xxz as xxz
import algorithm_lanczos as lanczos

# Parameters for Heisenberg model
L=12     # length of chain, keep it smaller than ~16, :-)
J=1      # strength of Heisenberg interaction
hs=0     # strangth of staggered magnetic field
sz=0     # magnetization, (number of up-spins = L//2 + sz)

rows, cols, data = stag.get_hamiltonian_sparse(L, J, hs, sz)
hamiltonian = sp.sparse.csr_matrix((data, (rows, cols)))
eigs, alphas, betas = lanczos.lanczos(hamiltonian)

print(eigs, alphas, betas)
