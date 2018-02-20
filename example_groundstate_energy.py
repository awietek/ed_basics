# -*- coding: utf-8 -*-
'''
Example script how to use the functions get_hamiltonian_sparse(...)
in the files:

- hamiltonian_tfi.py
- hamiltonian_hb_staggered.py
- hamiltonian_hb_xxz.py

:author: Alexander Wietek
:email: alexander.wietek@uibk.ac.at
:year: 2018
'''
from __future__ import absolute_import, division, print_function

import numpy as np
import scipy as sp
from scipy import sparse
from scipy import linalg
from scipy.sparse import linalg

import hamiltonian_tfi as tfi
import hamiltonian_hb_staggered as stag
import hamiltonian_hb_xxz as xxz

# Parameters for diagonalization
full_diagonalization = False
n_lowest_eigenvalues = 1

'''
Computation of the ground state energy of the transverse field Ising
model
'''
# Parameters for TFI model
L=12     # length of chain, keep it smaller than ~16, :-)
J=1      # strength of Ising interaction
hx=0.5   # strangth of transverse field

rows, cols, data = tfi.get_hamiltonian_sparse(L, J, hx)
hamiltonian = sp.sparse.csr_matrix((data, (rows, cols)))
if full_diagonalization:
    eigs = sp.linalg.eigh(hamiltonian.todense(), eigvals_only=True)
else:
    eigs = sp.sparse.linalg.eigsh(hamiltonian, k=n_lowest_eigenvalues,
                                  which='SA', return_eigenvectors=False,
                                  maxiter=1000)

print("Ground state energy of TFI model at J =", J, ", hx =",
      hx, ":", min(eigs))

'''
Computation of the ground state energy of the Heisenberg model
in a staggered magnetic field
'''
# Parameters for Heisenberg model
L=12     # length of chain, keep it smaller than ~16, :-)
J=1      # strength of Heisenberg interaction
hs=0.5   # strangth of staggered magnetic field
sz=0     # magnetization, (number of up-spins = L//2 + sz)

rows, cols, data = stag.get_hamiltonian_sparse(L, J, hs, sz)
hamiltonian = sp.sparse.csr_matrix((data, (rows, cols)))
if full_diagonalization:
    eigs = sp.linalg.eigh(hamiltonian.todense(), eigvals_only=True)
else:
    eigs = sp.sparse.linalg.eigsh(hamiltonian, k=n_lowest_eigenvalues,
                                  which='SA', return_eigenvectors=False,
                                  maxiter=1000)

print("Ground state energy of Heisenberg model with staggered field at J =",
      J, ", hs =", hs, ", Sz =", sz, ":", min(eigs))


'''
Computation of the ground state energy of the Heisenberg XXZ model
'''
# Parameters for Heisenberg model
L=12       # length of chain, keep it smaller than ~16, :-)
J=1        # strength of Heisenberg interaction
delta=0    # strength of XY anisotropy
sz=0       # magnetization, (number of up-spins = L//2 + sz)
k=1        # momentum quantum number
           #(in integer units, acutal momentum 2*pi/L*k

rows, cols, data = xxz.get_hamiltonian_sparse(L, J, delta, sz, k)
hamiltonian = sp.sparse.csr_matrix((data, (rows, cols)))
if full_diagonalization:
    eigs = sp.linalg.eigh(hamiltonian.todense(), eigvals_only=True)
else:
    eigs = sp.sparse.linalg.eigsh(hamiltonian, k=n_lowest_eigenvalues,
                                  which='SA', return_eigenvectors=False,
                                  maxiter=1000)
print("Ground state energy of Heisenberg model with staggered field at J =",
      J, ", delta =", delta, ", Sz =", sz, ", k = ", k, ":",  min(eigs))
