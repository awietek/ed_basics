# -*- coding: utf-8 -*-
'''
Example how to use the Lanczos algorithm in combination with the
sparse matrix data from creating Hamiltonians

:author: Alexander Wietek
:email: alexander.wietek@uibk.ac.at
:year: 2018
'''
from __future__ import absolute_import, division, print_function

import scipy as sp
import scipy.sparse

import hamiltonian_hb_xxz as xxz
import algorithm_lanczos as lanczos

# Parameters for Heisenberg model
L=12       # length of chain, keep it smaller than ~16, :-)
J=1        # strength of Heisenberg interaction
delta=0    # strength of XY anisotropy
sz=0       # magnetization, (number of up-spins = L//2 + sz)
k=1        # momentum quantum number
           #(in integer units, acutal momentum 2*pi/L*k

rows, cols, data = xxz.get_hamiltonian_sparse(L, J, delta, sz, k)
hamiltonian = sp.sparse.csr_matrix((data, (rows, cols)))
eigs, alphas, betas = lanczos.lanczos(hamiltonian)

print("Lanczos eigenvalues:", eigs)
print("alphas:", alphas)
print("betas:", betas)
