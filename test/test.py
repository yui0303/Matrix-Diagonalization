import sys
sys.path.append("../python")

import _matrix
import numpy as np
import time
from utils import *


def IsEqual(ANS, M, eps=1e-6):
    try:
        return np.allclose(np.array(ANS), np.array(M), atol=eps)
    except:
        return False

def is_valid(A, Q, R):
    if type(Q) != _matrix.Matrix: Q = _matrix.Matrix(Q.shape[0], Q.shape[1], Q.flatten())
    if type(R) != _matrix.Matrix: R = _matrix.Matrix(R.shape[0], R.shape[1], R.flatten())
    assert _matrix.is_orthogonal_matrix(Q)
    assert IsEqual(A, Q*R)
    
def create_matrix(nrow, ncol) -> (np.ndarray, _matrix.Matrix):
    # create symmetric matrix
    np_A = np.random.rand(nrow, ncol)
    # np_A = 0.5 * (np_A + np_A.T)
    np_list = np_A.flatten()
    A = _matrix.Matrix(nrow, ncol, np_list)
    return np_A, A

def test_constructor():
    np_A = np.random.rand(3, 3)
    np_list = np_A.flatten()
    A = _matrix.Matrix(3, 3, np_list)
    for i in range(3):
        for j in range(3):
            assert(A[i, j] == np_A[i][j])
    tmp_np_A = np.array(A)
    assert np.allclose(tmp_np_A, np_A)

def test_operator():
    np_A, A = create_matrix(3, 3)
    np_B, B = create_matrix(3, 3)
    assert np.allclose(A + B, np_A + np_B)
    assert np.allclose(A - B, np_A - np_B)
    assert np.allclose(A * B, np_A @ np_B)

def test_GS_QR_advance():
    for i in range(200):
        _, A = create_matrix(i+1, i+1)
        Q, R = _matrix.QR_decomposition(A)
        is_valid(A, Q, R)

def test_HS_QR_advance():
    for i in range(100):
        _, A = create_matrix(i+1, i+1)
        Q, R = _matrix.QR_decomposition(A, 2)
        is_valid(A, Q, R)

# def test_eigenvalue():
#     for i in range(100):
#         np_A, A = create_matrix(3, 3)
#         np_eig, _ = np.linalg.eig(np_A)
#         eig = _matrix.find_eigenvalue(A)
#         try:
#             assert np.allclose(sorted(np_eig), sorted(eig), atol=1e-3)
#         except:
#             print(A)
#             print(np_eig)
#             print(eig)

def test_base_qr():
    np_A = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    np_A = 0.5 * (np_A + np_A.T)
    np_list = np_A.flatten()
    A = _matrix.Matrix(3, 3, np_list)
    print()
    print('-'*80, end="\n")
    print("Test Base QR Decomposition A = QR")
    def wrapper(func, A, func_type):
        start = time.time()
        if func_type == 1: 
            Q, R = func(A, 2)
        else:
            Q, R = func(A)
        end = time.time()
        is_valid(A, Q, R)
        func_name = ["C++ GS", "C++ HS", "Numpy", "Custom Numpy GS", "Custom Numpy HS"]
        print("="*15 + " " + func_name[func_type] + " " + "="*15, end="\n")
        np_a = np.array(A)
        np_q = np.array(Q)
        np_r = np.array(R)
        for i in range(3):
            for j in range(3):
                print("{:10.5f}".format(np_a[i][j]), end="")
            print("\t", end="")
            for j in range(3):
                print("{:10.5f}".format(np_q[i][j]), end="")
            print("\t", end="")
            for j in range(3):
                print("{:10.5f}".format(np_r[i][j]), end="")
            print()
        print(f"Time cost: {(end - start)*1000:.4f} ms.")
    wrapper(_matrix.QR_decomposition, A, 0)
    wrapper(_matrix.QR_decomposition, A, 1)
    wrapper(np.linalg.qr, np_A, 2)
    wrapper(custom_QR_decomposition_GS, np_A, 3)
    wrapper(custom_QR_decomposition_HS, np_A, 4)
    print('-'*80, end="\n")

def benchmark():
    print("QR Decomposition Benchmark")
    print("{:>15s} {:>15s} {:>15s} {:>15s}".format("Size", "C++ GS", "C++ HS", "Numpy"))
    def wrapper(size, np_A, A):
        start = time.time()
        _, _ = _matrix.QR_decomposition(A)
        end = time.time()
        print("{:15d} {:15.4f} ".format(size, (end - start)*1000, 0), end="")
        start = time.time()
        if size <= 200:
            _, _ = _matrix.QR_decomposition(A, 2)
        end = time.time()
        print("{:15.4f}".format((end - start)*1000), end="")
        start = time.time()
        _, _ = np.linalg.qr(np_A)
        end = time.time()
        print("{:15.4f}".format((end - start)*1000), end="")
        print()

    size_list = [10, 15, 20, 30, 50, 100, 200, 500, 1000]
    for size in size_list:
        np_A, A = create_matrix(size, size)
        wrapper(size, np_A, A)
        
if __name__ == "__main__":
    test_base_qr()
    benchmark()