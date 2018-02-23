# -*- coding: utf-8 -*-
'''
Solution to Exercise 1a

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

# Step 1: just copied this function from hamiltonian_tfi.py
def get_hamiltonian_sparse(L, J, hx):
    '''
    Creates the Hamiltonian of the Transverse Field Ising model
    '''

    def get_site_value(state, site):
        ''' Function to get local value at a given site '''
        return (state >> site) & 1

    def hilbertspace_dimension(L):
        ''' return dimension of hilbertspace '''
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
            hamiltonian_rows.append(new_state)
            hamiltonian_cols.append(state)
            hamiltonian_data.append(hx)

    return hamiltonian_rows, hamiltonian_cols, hamiltonian_data


# Step 2: Compute the energies of the lowest 10 eigenvalues for x=0...2
# and plot results
L=12
J=1
n_lowest_eigenvalues = 10
for hx in np.linspace(0,2,20):
    rows, cols, data = get_hamiltonian_sparse(L, J, hx)
    hamiltonian = sp.sparse.csr_matrix((data, (rows, cols)))
    eigs = sp.sparse.linalg.eigsh(hamiltonian, k=n_lowest_eigenvalues,
                                  which='SA', return_eigenvectors=False,
                                  maxiter=1000)
    plt.plot(hx*np.ones_like(eigs), eigs-min(eigs), 's', color="steelblue")
    print(hx, eigs)
plt.xlabel(r"$h_x$")
plt.ylabel(r"$E - E_0$")
plt.ylim(-0.1, 6)
plt.legend()
plt.show()
