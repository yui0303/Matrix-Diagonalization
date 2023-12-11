#include "_eigen.hpp"

std::pair<Matrix, Matrix> QR_decomposition_wrapper(Matrix const& mat, int type)
{
    if (type == 1) {
        return QR_decomposition_GS(mat);
    } else if (type == 2) {
        return QR_decomposition_HS(mat);
    } else {
        throw std::invalid_argument("QR_decomposition_wrapper: invalid type");
    }
}

std::pair<Matrix, Matrix> QR_decomposition_GS(Matrix const& mat)
{
    Matrix Q = gram_schmidt(mat);
    Matrix R = Q.transpose() * mat;
    return {Q, R};
}

std::pair<Matrix, Matrix> QR_decomposition_HS(Matrix const& mat)
{
    Matrix Q = Matrix::Identity(mat.nrow(), mat.ncol());
    Matrix R = mat;
    Matrix H;
    for (size_t i = 0; i<mat.ncol()-1; ++i) 
    {
        std::vector<double> v(mat.nrow());
        for (size_t j = 0; j<mat.nrow(); ++j) 
        {
            if (j >= i)
                v[j] = R(j, i);
            else
                v[j] = 0;
        }
        H = householder(v, mat.nrow(), i);
        R = H * R;
        Q = Q * H;
    }
    return {Q, R};    
}

std::vector<double> find_eigenvalue(Matrix const& mat)
{
    size_t iteration = 5000000;
    
    Matrix A_old = mat;
    Matrix A_new = mat;
    double diff = std::numeric_limits<double>::max();
    double tolerance = 1e-10;
    size_t i = 0;
    while (diff > tolerance && i < iteration) {
        A_old = A_new;
        std::pair<Matrix, Matrix> QR = QR_decomposition_GS(A_old);
        A_new = QR.second * QR.first;
        diff = 0;
        for (size_t i = 0; i < A_old.nrow(); i++) {
            diff += std::abs(A_old(i, i) - A_new(i, i));
        }
        i++;
    }
    std::vector<double> eigenvalues;
    for (size_t i = 0; i < A_new.nrow(); i++) {
        eigenvalues.push_back(A_new(i, i)); // The eigenvalues are the diagonal elements
    }
    return eigenvalues;
}

std::vector<double> find_eigenvector(Matrix const& mat, double eigenvalue)
{
    return {};
}

bool is_orthogonal_matrix(Matrix const& mat)
{
    Matrix T = mat.transpose();
    Matrix I = mat * T;
    for (size_t i = 0; i < I.nrow(); i++) {
        for (size_t j = 0; j < I.ncol(); j++) {
            if (i == j) {
                if (std::abs(I(i, j) - 1) > EPSILON) {
                    return false;
                }
            } else {
                if (std::abs(I(i, j)) > EPSILON) {
                    return false;
                }
            }
        }
    }
    return true;
}
