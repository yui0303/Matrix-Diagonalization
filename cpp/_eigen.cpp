#include "_eigen.hpp"

std::pair<Matrix, Matrix> QR_decomposition_GS(Matrix const& mat)
{
    Matrix Q = gram_schmidt(mat);
    Matrix R = Q.transpose() * mat;
    return {Q, R};
}

std::pair<Matrix, Matrix> QR_decomposition_HS(Matrix const& mat)
{
    Matrix Q = householder(mat);
    Matrix R = Q.transpose() * mat;
    return {Q, R};
}

std::vector<double> find_eigenvalue(Matrix const& mat)
{
    std::pair<Matrix, Matrix> QR = QR_decomposition_GS(mat);
    Matrix A = QR.second * QR.first; // R*Q gives the next matrix in the sequence
    std::vector<double> eigenvalues;
    for (size_t i = 0; i < A.nrow(); i++) {
        eigenvalues.push_back(A(i, i)); // The eigenvalues are the diagonal elements
    }
    return eigenvalues;
}

std::vector<double> find_eigenvector(Matrix const& mat, double eigenvalue)
{
    return {};
}
