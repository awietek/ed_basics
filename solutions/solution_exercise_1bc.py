# -*- coding: utf-8 -*-
'''
Solution to Exercise 1b and 1c

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
import matplotlib.pyplot as plt

def get_site_value(state, site):
    ''' Function to get local value at a given site '''
    return (state >> site) & 1

def hilbertspace_dimension(L):
    ''' return dimension of hilbertspace '''
    return 2**L

def get_hamiltonian_sparse(L, J):
    ''' 
    Create Hamiltonian of Heisenberg model
    '''
    # Define chain lattice
    bonds = [(site, (site+1)%L) for site in range(L)]

    # Empty lists for sparse matrix
    hamiltonian_rows = []
    hamiltonian_cols = []
    hamiltonian_data = []

    # Run through all spin configurations
    for state in range(hilbertspace_dimension(L)):
        # Apply Ising bonds
        ising_diagonal = 0
        for bond in bonds:
            if get_site_value(state, bond[0]) == get_site_value(state, bond[1]):
                ising_diagonal += J/4
            else:
                ising_diagonal -= J/4
        hamiltonian_rows.append(state)
        hamiltonian_cols.append(state)
        hamiltonian_data.append(ising_diagonal)

        # Apply exchange bond
        for bond in bonds:
            flipmask = (1 << bond[0]) | (1 << bond[1])
            if get_site_value(state, bond[0]) != get_site_value(state, bond[1]):
                new_state = state ^ flipmask
                hamiltonian_rows.append(new_state)
                hamiltonian_cols.append(state)
                hamiltonian_data.append(0.5*J)
                
    return hamiltonian_rows, hamiltonian_cols, hamiltonian_data

def get_spinspin_sparse(L, r, sz_only=True):
    ''' 
    Create operator S_0 S_r (if sz_only=False) or S^z_0 S^z_r
    if sz_only = True
    '''
    # Empty lists for sparse matrix
    rows = []
    cols = []
    data = []

    bond = (0,r)
    
    # Run through all spin configurations
    for state in range(hilbertspace_dimension(L)):
        # Apply Ising bonds
        ising_diagonal = 0
        if get_site_value(state, bond[0]) == get_site_value(state, bond[1]):
            ising_diagonal += J/4
        else:
            ising_diagonal -= J/4
        rows.append(state)
        cols.append(state)
        data.append(ising_diagonal)

        if not sz_only:
            # Apply exchange bond
            flipmask = (1 << bond[0]) | (1 << bond[1])
            if get_site_value(state, bond[0]) != get_site_value(state, bond[1]):
                new_state = state ^ flipmask
                rows.append(new_state)
                cols.append(state)
                data.append(0.5)

    return rows, cols, data

# 12 site antiferro Heisenberg chain
L=12
J=1
print("system size L =", L)

# Get the ground state
h_rows, h_cols, h_data = get_hamiltonian_sparse(L, J)
hamiltonian = sp.sparse.csr_matrix((h_data, (h_rows, h_cols)))
energy, groundstate = sp.sparse.linalg.eigsh(hamiltonian, k=1,
                                             which='SA',
                                             return_eigenvectors=True,
                                             maxiter=1000)
print("ground state energy (per site):", energy/L)

# Compute the ground state
groundstate = groundstate[:,0]
sz_corrs = []
spin_corrs = []
for r in range(1,L,1):
    sz_rows, sz_cols, sz_data = get_spinspin_sparse(L,r,True)
    sz_matrix = sp.sparse.csr_matrix((sz_data, (sz_rows, sz_cols)))
    sz_corrs.append(np.dot(groundstate, sz_matrix.dot(groundstate)))
    spin_rows, spin_cols, spin_data = get_spinspin_sparse(L,r,False)
    spin_matrix = sp.sparse.csr_matrix((spin_data, (spin_rows, spin_cols)))
    spin_corrs.append(np.dot(groundstate, spin_matrix.dot(groundstate)))

# Check 3* <S^z_0S^z_r> = <S_0 S_r> 
assert np.allclose(3*np.array(sz_corrs), np.array(spin_corrs))

# Plot the results
plt.plot(range(1,L,1), sz_corrs, "s-", label=r"$C(r)=\langle 0 | S^z_0 S^z_r | 0 \rangle$")
plt.plot(range(1,L,1), spin_corrs, "s-", label=r"$C(r)=\langle 0 | \mathbf{S}_0\cdot \mathbf{S}_r | 0 \rangle$")
plt.xlabel(r"$r$")
plt.ylabel(r"$C(r)$")
plt.legend()
plt.show()
