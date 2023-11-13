#include "_operation.hpp"

double dot_product(const std::vector<double>& a, const std::vector<double>& b) 
{
    double result = 0;
    for (size_t i = 0; i < a.size(); i++) {
        result += a[i] * b[i];
    }
    return result;
}

std::vector<double> subtract_vector(const std::vector<double> &a, const std::vector<double> &b)
{
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        result[i] = a[i] - b[i];
    }
    return result;
}

std::vector<double> scalar_vector(const std::vector<double> &a, double scalar)
{
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        result[i] = a[i] * scalar;
    }
    return result;
}

std::vector<double> subtract_projection(const std::vector<double>& a, const std::vector<double>& b) 
{
    double projection = dot_product(a, b) / dot_product(b, b);
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        result[i] = a[i] - projection * b[i];
    }
    return result;
}

std::vector<double> normalize(const std::vector<double>& a) 
{
    double norm = sqrt(dot_product(a, a));
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        result[i] = a[i] / norm;
    }
    return result;
}

Matrix gram_schmidt(Matrix const& mat) 
{
    Matrix Q(mat.nrow(), mat.ncol());
    Matrix TansMat = mat.transpose();
    std::cout<< TansMat<<std::endl;
    for (size_t i = 0; i < mat.ncol(); i++) {
        std::vector<double> vi = TansMat(i);
        // std::cout<< TansMat(i)<<std::endl;
        for(size_t j = 0; j < i; j++) {
            //std::cout<<Q(j)<<std::endl;
            double scalar = dot_product(TansMat(i), Q(j))/dot_product(Q(j), Q(j));
            vi = subtract_vector(vi, scalar_vector(Q(j), scalar));
        }
        vi = normalize(vi);
        for(size_t j = 0; j < vi.size(); j++) {
            Q(i, j) = vi[j];
        }
    }
    return Q.transpose();
}

//function to multiply two matrix m*n & n*p with tile size and m, n, p are all multiples of tile size
Matrix multiply_tile(Matrix const& mat1, Matrix const& mat2, size_t tsize)
{
    //Check whether the multiplication is possible
    if (mat1.ncol() != mat2.nrow())
    {
        throw std::out_of_range("invalid matrix dimensions for multiplication");
    }

    // Create a new matrix to store the result of the multiplication.
    Matrix result(mat1.nrow(), mat2.ncol());
    size_t m = mat1.nrow();
    size_t n = mat2.ncol();
    size_t p = mat1.ncol();
    // Divide the matrices into tiles of size tile_m x tile_n and tile_n x tile_p.
    for (size_t i = 0; i < m; i += tsize)
    {
        for (size_t j = 0; j < n; j += tsize)
        {
            for (size_t k = 0; k < p; k += tsize)
            {
                // Multiply each tile of mat1 with the corresponding tile of mat2 and add the result to the corresponding tile of the result matrix.
                for (size_t ii = i; ii < std::min(i + tsize, m); ++ii)
                {
                    for (size_t jj = j; jj < std::min(j + tsize, n); ++jj) 
		    {
		    	double sum = .0;
                        for (size_t kk = k; kk < std::min(k + tsize, p); ++kk)
                        {
                            sum += mat1(ii, kk) * mat2(kk, jj);
                        }
			result(ii, jj) += sum;
                    }
                }
            }
        }
    }

    return result;
}

Matrix multiply_naive(Matrix const& mat1, Matrix const& mat2)
{
    // Check if the dimensions of the matrices are valid for multiplication.
    if (mat1.ncol() != mat2.nrow())
    {
        throw std::out_of_range("invalid matrix dimensions for multiplication");
    }

    // Create a new matrix to store the result of the multiplication.
    Matrix result(mat1.nrow(), mat2.ncol());

    // Multiply each element of mat1 with the corresponding element of mat2 and add the result to the corresponding element of the result matrix.
    for (size_t i = 0; i < mat1.nrow(); ++i)
    {
        for (size_t j = 0; j < mat2.ncol(); ++j)
        {
	    double sum = .0;
            for (size_t k = 0; k < mat1.ncol(); ++k)
            {
                sum += mat1(i, k) * mat2(k, j);
            }
	    result(i, j) = sum;
        }
    }

    return result;
}


/*
Matrix multiply_mkl(Matrix const& mat1, Matrix const& mat2)
{
    // Check if the dimensions of the matrices are valid for multiplication.
    if (mat1.ncol() != mat2.nrow())
    {
        throw std::out_of_range("invalid matrix dimensions for multiplication");
    }

    // Create a new matrix to store the result of the multiplication.
    Matrix result(mat1.nrow(), mat2.ncol());

    // Compute the matrix multiplication using DGEMM.
    cblas_dgemm(
        CblasRowMajor, 
        CblasNoTrans, 
        CblasNoTrans, 
        mat1.nrow(), 
        mat2.ncol(), 
        mat1.ncol(), 
        1.0, 
        mat1.data(), 
        mat1.ncol(), 
        mat2.data(), 
        mat2.ncol(), 
        0.0, 
        result.data(), 
        result.ncol()
        );

    return result;
}
*/