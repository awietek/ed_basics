# -*- coding: utf-8 -*-
'''
Solution to Exercise 2b

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


# Step 1: change Hamiltonian of TFI model to work with spinflip symmetry
# changes are marked in CAPITALS and #s
def get_hamiltonian_sparse(L, J, hx, gamma):  # INTRODUCE GAMMA AS PARAMETER #
    '''
    Creates the Hamiltonian of the Transverse Field Ising model

    gamma:  +1 or -1 depending if even or odd representation of Z2 spinflip
           symmetry is chosen
    '''

    def get_site_value(state, site):
        ''' Function to get local value at a given site '''
        return (state >> site) & 1

    def hilbertspace_dimension(L):
        ''' return dimension of hilbertspace '''
        return 2**L

    # INTRODUCED FLIP FUNCTION                   #
    def flip(L, state):                          #
        ''' flip all spins in a given state '''  #
        return ((1 << L)-1) & (~state)           #

    # INTRODUCED GET_REPRESENTATIVE                                      #
    def get_representative(L, state):                                    #
        ''' finds representative and representative flip for a state ''' #
        flipped_state = flip(L, state)                                   #
        if flipped_state < state:                                        #
            return flipped_state, gamma                                   #
        else:                                                            #
            return state, 1                                              #

    # CREATE REPRESENTATIVES, (NO NORMS NEEDED FOR SPINFLIP SYMMETRY)    #
    basis_states = []                                                    # 
    for state in range(hilbertspace_dimension(L)):                       # 
        representative, factor = get_representative(L, state)            #
        if state == representative:                                      #
            basis_states.append(state)                                   #
            
    
    # Define chain lattice
    ising_bonds = [(site, (site+1)%L) for site in range(L)]

    # Empty lists for sparse matrix
    hamiltonian_rows = []
    hamiltonian_cols = []
    hamiltonian_data = []

    # Run through all spin configurations
    # for state in range(hilbertspace_dimension(L)):
    for state_index, state in enumerate(basis_states): # REPLACED PREVIOUS LINE #
        # Apply Ising bonds
        ising_diagonal = 0
        for bond in ising_bonds:
            if get_site_value(state, bond[0]) == get_site_value(state, bond[1]):
                ising_diagonal += J
            else:
                ising_diagonal -= J
        hamiltonian_rows.append(state)
        hamiltonian_cols.append(state)
        hamiltonian_data.append(ising_diagonal)

        # Apply transverse field
        for site in range(L):
            # Flip spin at site
            new_state = state ^ (1 << site)
            # GET THE REPRESENTATIVE AND ITS INDEX                    #
            representative, factor = get_representative(L, new_state) #
            new_state_index = basis_states.index(representative)      #

            # USE INDICES INSTEAD OF STATED HERE       #
            hamiltonian_rows.append(new_state_index)   #          
            hamiltonian_cols.append(state_index)       #
            hamiltonian_data.append(hx*factor)         #

    return hamiltonian_rows, hamiltonian_cols, hamiltonian_data

# Step 2: Compute the energies of the lowest 10 eigenvalues for x=0...2
# and bot symmetry sectors gamma=1,-1 and plot results
L=12
J=1
n_lowest_eigenvalues = 10
for hx in np.linspace(0,2,20):
    eigs = []
    for gamma in [1, -1]:
        rows, cols, data = get_hamiltonian_sparse(L, J, hx, gamma)
        hamiltonian = sp.sparse.csr_matrix((data, (rows, cols)))
        eigs.append(sp.sparse.linalg.eigsh(hamiltonian,
                                           k=n_lowest_eigenvalues,
                                           which='SA',
                                           return_eigenvectors=False,
                                           maxiter=1000))
    groundstate_energy = min(min(eigs[0]),min(eigs[1]))
    if hx==0:
        plt.plot(hx*np.ones_like(eigs[0]), eigs[0]-groundstate_energy,
                 's', color="steelblue", label="even spinflip")
        plt.plot(hx*np.ones_like(eigs[1]), eigs[1]-groundstate_energy,
                 's', color="firebrick", label="odd spinflip")
    else:
        plt.plot(hx*np.ones_like(eigs[0]), eigs[0]-groundstate_energy,
                 's', color="steelblue")
        plt.plot(hx*np.ones_like(eigs[1]), eigs[1]-groundstate_energy,
                 's', color="firebrick")
    print(hx, eigs)
plt.xlabel(r"$h_x$")
plt.ylabel(r"$E - E_0$")
plt.ylim(-0.1, 6)
plt.legend()
plt.show()
