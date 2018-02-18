# -*- coding: utf-8 -*-
"""
Basic script for Exact Diagonalization of transverse field Ising model

:author: Alexander Wietek
"""
from __future__ import division, print_function
import numpy as np
import scipy as sp
from scipy import sparse
from scipy import linalg
from scipy.sparse import linalg

# Parameters for TFI model
L=12     # length of chain, keep it smaller than ~16, :-)
J=1      # strength of Ising interaction
hx=1     # strangth of transverse field

# Parameters for diagonalization
full_diagonalization = False
n_lowest_eigenvalues = 5

def get_site_value(state, site):
    ''' Function to get local value at a given site '''
    return (state >> site) & 1 

def hilbertspace_dimension(L):
    return 2**L

def get_hamiltonian_sparse(L, J, hx):
    # Define chain lattice
    ising_bonds = [(site, (site+1)%L) for site in range(L)]

    # Empty lists for sparse matrix
    hamiltonian_rows = []
    hamiltonian_cols = []
    hamiltonian_data = []
    
    # Run through all spin configurations
    for state in range(hilbertspace_dimension(L)):
        
        # Apply Ising bonds
        ising_diagonal = 0
        for bond in ising_bonds:
            if get_site_value(state, bond[0]) == get_site_value(state, bond[1]):
                ising_diagonal += J/4
            else:
                ising_diagonal -= J/4
        hamiltonian_rows.append(state)
        hamiltonian_cols.append(state)
        hamiltonian_data.append(ising_diagonal)
            
        # Apply transverse field
        for site in range(L):

            # Flip spin at site 
            new_state = state ^ (1 << site)
            hamiltonian_rows.append(new_state)
            hamiltonian_cols.append(state)
            hamiltonian_data.append(hx / 2)

            hamiltonian_rows.append(state)
            hamiltonian_cols.append(new_state)
            hamiltonian_data.append(hx / 2)
    return hamiltonian_rows, hamiltonian_cols, hamiltonian_data

rows, cols, data = get_hamiltonian_sparse(L, J, hx)
hamiltonian = sp.sparse.csr_matrix((data, (rows, cols)))
if hilbertspace_dimension(L) <= n_lowest_eigenvalues or full_diagonalization:
    eigs = sp.linalg.eigh(hamiltonian.todense(), eigvals_only=True)
else:
    eigs = sp.sparse.linalg.eigsh(hamiltonian, k=n_lowest_eigenvalues,
                                  which='SA', return_eigenvectors=False,
                                  maxiter=1000)
print(eigs)
