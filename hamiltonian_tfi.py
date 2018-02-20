# -*- coding: utf-8 -*-
"""
Function to create the Hamiltonian of the transverse field Ising model,
no symmetries implemented

:author: Alexander Wietek
:email: alexander.wietek@uibk.ac.at
:year: 2018
"""
from __future__ import absolute_import, division, print_function

import numpy as np
import scipy as sp
from scipy import sparse
from scipy import linalg
from scipy.sparse import linalg

def get_hamiltonian_sparse(L, J, hx):
    '''
    Creates the Hamiltonian of the Transverse Field Ising model
    on a linear chain lattice with periodic boundary conditions.
       
    Args:
        L (int): length of chain
        J (float): coupling constant for Ising term
        hx (float): coupling constant for transverse field
    
    Returns:
        (hamiltonian_rows, hamiltonian_cols, hamiltonian_data) where:
        hamiltonian_rows (list of ints): row index of non-zero elements
        hamiltonian_cols (list of ints): column index of non-zero elements
        hamiltonian_data (list of floats): value of non-zero elements
    '''
    
    def get_site_value(state, site):
        ''' Function to get local value at a given site '''
        return (state >> site) & 1 

    def hilbertspace_dimension(L):
        return 2**L

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

    return hamiltonian_rows, hamiltonian_cols, hamiltonian_data
