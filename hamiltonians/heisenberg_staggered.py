# -*- coding: utf-8 -*-
"""
Basic script for Exact Diagonalization of spin 1/2 Heisenberg model
with staggered magnetic field

:author: Alexander Wietek
"""
from __future__ import division, print_function
import numpy as np
import scipy as sp
from scipy import sparse
from scipy import linalg
from scipy.sparse import linalg

# Parameters for Heisenberg model
L=12     # length of chain, keep it smaller than ~16, :-)
J=1      # strength of Heisenberg interaction
hs=0     # strangth of staggered magnetic field
sz=0     # magnetization, (number of up-spins = L//2 + sz)

# Parameters for diagonalization
full_diagonalization = False
n_lowest_eigenvalues = 5

# Functions to manipulate states
def get_site_value(state, site):
    ''' Function to get local value at a given site '''
    return (state >> site) & 1 

def set_site_value(state, site, value):
    ''' Function to set local value at a given site '''
    return (state & ~(1 << site)) | (value << site)


# Functions to create sz basis
def first_state(L, sz):
    ''' Return first state of Hilbert space in lexicographic order '''
    n_upspins = L//2 + sz
    return (1 << n_upspins) - 1

def next_state(state):
    ''' 
    Return next state of Hilbert space in lexicographic order 

    This function implements is a nice trick for spin 1/2 only,
    see http://graphics.stanford.edu/~seander/bithacks.html
    #NextBitPermutation for details
    '''
    t = (state | (state - 1)) + 1;  
    return t | ((((t & -t) // (state & -state)) >> 1) - 1)

def last_state(L, sz):
    ''' Return last state of Hilbert space in lexicographic order '''
    n_upspins = L//2 + sz
    return ((1 << n_upspins) - 1) << (L - n_upspins)

# Function to create Hamiltonian
def get_hamiltonian_sparse(L, J, hs, sz):
    # check if sz is valid
    assert((sz <= (L // 2 + L % 2)) and (sz >= -L//2))

    # Create list of states with fixed sz
    basis_states = []
    state = first_state(L, sz)
    end_state = last_state(L, sz)
    while state <= end_state:
        basis_states.append(state)
        state = next_state(state)

    # Define chain lattice
    heisenberg_bonds = [(site, (site+1)%L) for site in range(L)]

    # Empty lists for sparse matrix
    hamiltonian_rows = []
    hamiltonian_cols = []
    hamiltonian_data = []

    # Run through all spin configurations
    for state_index, state in enumerate(basis_states):

        # Apply diagonal Ising bonds
        diagonal = 0
        for bond in heisenberg_bonds:
            if get_site_value(state, bond[0]) == get_site_value(state, bond[1]):
                diagonal += J/4
            else:
                diagonal -= J/4
                
        # Apply diagonal staggered Sz field
        for site in range(0, L, 2):
            diagonal += hs*(2*get_site_value(state, site) - 1)
            diagonal -= hs*(2*get_site_value(state, site+1) - 1)

        hamiltonian_rows.append(state_index)
        hamiltonian_cols.append(state_index)
        hamiltonian_data.append(diagonal)
        
        # Apply exchange interaction
        for bond in heisenberg_bonds:
            flipmask = (1 << bond[0]) | (1 << bond[1])
            if get_site_value(state, bond[0]) != get_site_value(state, bond[1]):
                new_state = state ^ flipmask
                new_state_index = basis_states.index(new_state)
                hamiltonian_rows.append(state_index)
                hamiltonian_cols.append(new_state_index)
                hamiltonian_data.append(J/2)
                hamiltonian_rows.append(new_state_index)
                hamiltonian_cols.append(state_index)
                hamiltonian_data.append(J/2)

    return hamiltonian_rows, hamiltonian_cols, hamiltonian_data
    
rows, cols, data = get_hamiltonian_sparse(L, J, hs, sz)
hamiltonian = sp.sparse.csr_matrix((data, (rows, cols)))

eigs = sp.sparse.linalg.eigsh(hamiltonian, k=n_lowest_eigenvalues,
                              which='SA', return_eigenvectors=False,
                              maxiter=1000)
print(eigs)
