# -*- coding: utf-8 -*-
'''
Implementation of basic Lanczos algorithm

:author: Alexander Wietek
:email: alexander.wietek@uibk.ac.at
:year: 2018
'''
from __future__ import absolute_import, division, print_function

import numpy as np
import scipy as sp
from scipy import sparse

def lanczos(matrix, precision=1e-12, max_iterations = 1000):
    '''
    Basic Lanczos algorithm without reorthogonalization
       
    Args:
        matrix (scipy.sparse.csr_matrix): Matrix for Lanczos iterations
        precision (float): precision for convergence of minimal eigenvalue
        max_iterations (int): maximum amount of Lanczos steps
    
    Returns:
        (t_eigs, alphas, betas) where:
        t_eigs (numpy.array): final eigenvalues of T-matrix
        alphas (numpy.array): Diagonal elements of the T-matrix
        betas (numpy.array): Off-diagonal elements of the T-Matrix
    '''   
    linear_size = matrix.get_shape()[0]
    
    alphas = []
    betas = [0.]
    eigenvalues = []

    v1 = np.random.rand(linear_size) + 1j*np.random.rand(linear_size)
    v1 /= np.linalg.norm(v1)
    v0 = np.zeros(linear_size, dtype=complex)
    w = np.zeros(linear_size, dtype=complex)
    
    alpha = 0.
    beta = 0.

    prev_energy = 0
    for iteration in range(1, max_iterations):

        # Perform basic Lanczos iteration step
        w = matrix.dot(v1)
        alpha = np.real(np.dot(np.conj(v1),w))
        w = w - alpha*v1 - beta*v0
        v0 = np.copy(v1)
        beta = np.real(np.sqrt(np.dot(np.conj(w),w)))
        v1 = 1/beta*w
        alphas.append(alpha)
        betas.append(beta)

        # Build up T-matrix
        t_matrix = np.diag(np.array(alphas)) + \
                   np.diag(np.array(betas)[1:-1],k=1) + \
                   np.diag(np.array(betas)[1:-1],k=-1)
        t_eigs = np.linalg.eigvalsh(t_matrix)

        # Return if converged
        if np.abs(min(t_eigs) - prev_energy) < precision:
            print("Lanczos converged in", iteration, "steps")
            return t_eigs, np.array(alphas), np.array(betas[1:])
        prev_energy = min(t_eigs)
