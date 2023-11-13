#ifndef __EIGEN_HPP__
#define __EIGEN_HPP__

#include "_matrix.hpp"
#include "_operation.hpp"

std::pair<Matrix, Matrix> QR_decomposition(Matrix const& mat);
std::vector<double> find_eigenvalue(Matrix const& mat);
std::vector<double> find_eigenvector(Matrix const& mat, double eigenvalue);

#endif