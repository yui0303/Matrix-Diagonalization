'''
householder, custom_QR_decomposition_HS is taken from https://rosettacode.org/wiki/QR_decomposition#Python
'''

import numpy as np

def null_space(mat, valid_row):
    # Step 1: Create an augmented matrix [A|0]
    augmented_mat = np.zeros((valid_row, mat.shape[1] * 2))
    augmented_mat[:, :mat.shape[1]] = mat[:valid_row, :]

    # Step 2: Perform Gaussian elimination
    lead = 0
    for r in range(valid_row):
        if lead >= mat.shape[1]:
            return []
        i = r
        while augmented_mat[i, lead] == 0:
            i += 1
            if i == valid_row:
                i = r
                lead += 1
                if mat.shape[1] == lead:
                    return []
        # Swap rows i and r
        augmented_mat[[i, r]] = augmented_mat[[r, i]]
        lv = augmented_mat[r, lead]
        augmented_mat[r] /= lv
        for i in range(valid_row):
            if i != r:
                sub = augmented_mat[i, lead]
                augmented_mat[i] -= augmented_mat[r] * sub
        lead += 1

    # Step 3: Identify pivot columns
    isPivot = np.zeros(mat.shape[1], dtype=bool)
    for r in range(valid_row):
        for c in range(mat.shape[1]):
            if augmented_mat[r, c] == 1:
                isPivot[c] = True
                break

    # Step 4: Solve the system for each free variable
    null_space_vectors = []
    for i in range(len(isPivot)):
        if not isPivot[i]:
            special_solution = np.zeros(mat.shape[1])
            special_solution[i] = 1
            for r in range(valid_row):
                if augmented_mat[r, i] != 0:
                    special_solution[r] = -augmented_mat[r, i]
            null_space_vectors.append(special_solution)

    # Return the null space
    return null_space_vectors

def gram_schmidt(mat):
    Q = np.zeros_like(mat)
    TansMat = mat.T
    num_non_zero_vec = 0
    for i in range(mat.shape[1]):
        vi = TansMat[i]
        for j in range(num_non_zero_vec):
            scalar = np.dot(TansMat[i], Q[j]) / np.dot(Q[j], Q[j])
            vi = vi - scalar * Q[j]
        if not np.allclose(vi, 0):
            vi = vi / np.linalg.norm(vi)
            Q[num_non_zero_vec] = vi
            num_non_zero_vec += 1
    null_vecs = null_space(Q, num_non_zero_vec)
    
    for i in range(num_non_zero_vec, mat.shape[1]):
        vi = null_vecs[i - num_non_zero_vec]
        vi = vi / np.linalg.norm(vi)
        Q[i] = vi
    return Q.T

def householder(A):
    v = A / (A[0] + np.copysign(np.linalg.norm(A), A[0]))
    v[0] = 1
    H = np.eye(A.shape[0])
    H -= (2 / np.dot(v, v)) * np.dot(v[:, None], v[None, :])
    return H

def custom_QR_decomposition_GS(A):
    Q = gram_schmidt(A)
    R = Q.T @ A
    return Q, R

def custom_QR_decomposition_HS(A):
    n = A.shape[0]
    Q = np.eye(n)
    for i in range(n - 1):
        H = np.eye(n)
        H[i:, i:] = householder(A[i:, i])
        Q = Q @ H
        A = H @ A
    return Q, A
    