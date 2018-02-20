# -*- coding: utf-8 -*-
"""
Function to create the Hamiltonian of the Heisenberg XXZ model with Ising
anisotropy, Sz and momentum conservation implemented

:author: Alexander Wietek
:email: alexander.wietek@uibk.ac.at
:year: 2018
"""
from __future__ import division, print_function
import numpy as np

def get_hamiltonian_sparse(L, J, delta, sz, k):
    '''
    Creates the Hamiltonian of the Heisenberg XXZ model with Ising
    anisotropy on a linear chain lattice with periodic boundary conditions.

    Args:
        L (int): length of chain
        J (float): coupling constant for Heisenberg term
        delta (float): strength of Ising anisotropy
        sz (int): total Sz
        k (int): defines momentum sector (k = 0,...,L-1)
                 actual momentum is given by 2*pi*k / L
    Returns:
        (hamiltonian_rows, hamiltonian_cols, hamiltonian_data) where:
        hamiltonian_rows (list of ints): row index of non-zero elements
        hamiltonian_cols (list of ints): column index of non-zero elements
        hamiltonian_data (list of floats): value of non-zero elements
    '''

    # Functions to manipulate states
    def get_site_value(state, site):
        ''' Function to get local value at a given site '''
        return (state >> site) & 1

    def set_site_value(state, site, value):
        ''' Function to set local value at a given site '''
        site_val = (value << site)
        return (state ^ site_val) | site_val

    # Functions to create sz basis and compute translations
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
        t = (state | (state - 1)) + 1
        return t | ((((t & -t) // (state & -state)) >> 1) - 1)

    def last_state(L, sz):
        ''' Return last state of Hilbert space in lexicographic order '''
        n_upspins = L//2 + sz
        return ((1 << n_upspins) - 1) << (L - n_upspins)

    def translate(L, state, n_translation_sites):
        ''' translates a state by n_translation_sites '''
        new_state = 0
        for site in range(L):
            site_value = get_site_value(state, site)
            new_state = set_site_value(new_state, (site + n_translation_sites)%L, site_value)
        return new_state

    def get_representative(L, state):
        ''' finds representative and representative translation for a state '''
        representative = state
        translation = 0
        for n_translation_sites in range(L):
            new_state = translate(L, state, n_translation_sites)
            if new_state < representative:
                representative = new_state
                translation = n_translation_sites
        return representative, translation


    # check if sz is valid
    assert (sz <= (L // 2 + L % 2)) and (sz >= -L//2)

    # Create list of representatives and norms
    basis_states = []
    norms = []

    state = first_state(L, sz)
    end_state = last_state(L, sz)
    while state <= end_state:
        representative, translation = get_representative(L, state)
        if state == representative:

            # Compute the normalization constant
            amplitude = 0.0
            for n_translation_sites in range(L):
                new_state = translate(L, state, n_translation_sites)
                if new_state == state:
                    amplitude += np.exp(1j*2*np.pi*k/L*n_translation_sites)
            norm = np.sqrt(np.abs(amplitude))
            if norm > 1e-12:
                basis_states.append(state)
                norms.append(norm)

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
        hamiltonian_rows.append(state_index)
        hamiltonian_cols.append(state_index)
        hamiltonian_data.append(diagonal)

        # Apply exchange interaction
        for bond in heisenberg_bonds:
            flipmask = (1 << bond[0]) | (1 << bond[1])
            if get_site_value(state, bond[0]) != get_site_value(state, bond[1]):
                new_state = state ^ flipmask
                representative, translation = get_representative(L, new_state)
                new_state_index = basis_states.index(representative)
                coeff = ((J+delta)/2)*(norms[new_state_index]/\
                                       norms[state_index])*\
                                       np.exp(1j*2*np.pi*k/L*translation)/2
                hamiltonian_rows.append(new_state_index)
                hamiltonian_cols.append(state_index)
                hamiltonian_data.append(coeff)
                hamiltonian_rows.append(state_index)
                hamiltonian_cols.append(new_state_index)
                hamiltonian_data.append(np.conj(coeff))

    return hamiltonian_rows, hamiltonian_cols, hamiltonian_data
