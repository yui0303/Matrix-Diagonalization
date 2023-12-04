import sys
# sys.path.append('../cpp')
import _matrix
import numpy as np
import time

def IsEqual(ANS, M, eps=1e-6):
    if ANS.nrow != M.nrow or ANS.ncol != M.ncol:
        return False
    for i in range(ANS.nrow):
        for j in range(ANS.ncol):
            if abs(ANS[i, j] - M[i, j]) > eps:
                return False
    return True

def create_matrix(nrow, ncol) -> (np.ndarray, _matrix.Matrix):
    np_A = np.random.rand(nrow, ncol)
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

def test_base():
    np_A = np.random.rand(3, 3)
    np_A = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    #np_A = np.array([[1, 2, 0], [1, 2, 0], [1, 2, 0]])

    np_list = np_A.flatten()
    
    print("=========== Numpy ===========")
    print(np_A)
    np_Q, np_R = np.linalg.qr(np_A)
    print("numpy Q:")
    print(np_Q)
    print("numpy R:")
    print(np_R)
    print("========== Q*R ==========")
    B_np = np.matmul(np_Q, np_R)
    print(B_np)
    print("========== _matrix ==========")
    A = _matrix.Matrix(3, 3, np_list)
    print(A)
    Q, R = _matrix.QR_decomposition(A)
    print("_matrix Q:")
    print(Q)
    print("_matrix R:")
    print(R)
    print("========== Q*R ==========")
    B = _matrix.Matrix(3, 3)
    B = Q*R
    print(B)

    assert(IsEqual(A, B))

def test_advance():
    for i in range(100):
        _, A = create_matrix(i+1, i+1)
        Q, R = _matrix.QR_decomposition(A)
        B = Q*R
        assert(IsEqual(A, B))
def test_1000():
    _, A = create_matrix(1000, 1000)
    Q, R = _matrix.QR_decomposition(A)
    B = Q*R
    assert(IsEqual(A, B))

def run_wrapper(func, A, func_type, size):
    def wrapper():
        start = time.time()
        _, _ = func(A)
        end = time.time()
        print(f"{func_type} function took {end - start:.4f} seconds to run ({size}).")
    return wrapper

def benchmark():
    size_list = [100, 200, 500, 1000]
    for size in size_list:
        np_A, A = create_matrix(size, size)
        run_wrapper(_matrix.QR_decomposition, A, "C++", size)()
        run_wrapper(np.linalg.qr, np_A, "Numpy", size)()
if __name__ == "__main__":
    test_base()
    benchmark()



    