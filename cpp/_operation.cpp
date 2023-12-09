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

std::vector<double> scalar_vector(double scalar, const std::vector<double> &a)
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

bool is_zero_vector(const std::vector<double>& a) 
{
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i] > EPSILON) {
            return false;
        }
    }
    return true;
}

Matrix gram_schmidt(Matrix const& mat) 
{
    Matrix Q(mat.nrow(), mat.ncol());
    Matrix TansMat = mat.transpose();
    size_t num_non_zero_vec = 0;
    // bool sign = false;
    for (size_t i = 0; i < mat.ncol(); i++) {
        std::vector<double> vi = TansMat(i);
        for(size_t j = 0; j < num_non_zero_vec; j++) {
            double scalar = dot_product(TansMat(i), Q(j))/dot_product(Q(j), Q(j));
            vi = subtract_vector(vi, scalar_vector(scalar, Q(j)));
        }
        if(!is_zero_vector(vi)) {
            vi = normalize(vi);
            for(size_t j = 0; j < vi.size(); j++) {
                Q(num_non_zero_vec, j) = vi[j];
            }
            ++num_non_zero_vec;
        }
    }

    std::vector<std::vector<double>> null_vecs = null_space(Q, num_non_zero_vec);
    
    for (size_t i=num_non_zero_vec; i<mat.ncol(); i++) {
        std::vector<double> vi = null_vecs[i-num_non_zero_vec];
        vi = normalize(vi);
        for (size_t j=0; j<mat.nrow(); j++) {
            Q(i, j) = vi[j];
        }
    }
    return Q.transpose();
}

std::vector<std::vector<double>> null_space(Matrix const& mat, size_t valid_row)
{
    // Step 1: Create an augmented matrix [A|0]
    Matrix augmented_mat(valid_row, mat.ncol() * 2);
    for (size_t i = 0; i < valid_row; i++) {
        for (size_t j = 0; j < mat.ncol(); j++) {
            augmented_mat(i, j) = mat(i, j);
        }
    }
    // std::cout<<"augmented_mat: \n"<<augmented_mat<<std::endl;
    // std::cout<<"pass1"<<std::endl;
    // Step 2: Perform Gaussian elimination
    size_t lead = 0;
    for (size_t r = 0; r < valid_row; r++) {
        if (lead >= mat.ncol())
            return {};
        size_t i = r;
        while (augmented_mat(i, lead) == 0) {
            i++;
            if (i == valid_row) {
                i = r;
                lead++;
                if (mat.ncol() == lead)
                    return {};
            }
        }
        // Swap rows i and r
        for (size_t j = 0; j < augmented_mat.ncol(); j++) {
            std::swap(augmented_mat(i, j), augmented_mat(r, j));
        }
        auto lv = augmented_mat(r, lead);
        for (size_t j = 0; j < augmented_mat.ncol(); j++) {
            augmented_mat(r, j) /= lv;
        }
        for (size_t i = 0; i < valid_row; i++) {
            if (i != r) {
                auto sub = augmented_mat(i, lead);
                for (size_t j = 0; j < augmented_mat.ncol(); j++) {
                    augmented_mat(i, j) -= (augmented_mat(r, j) * sub);
                }
            }
        }
        lead++;
    }
    // std::cout<<"pass2"<<std::endl;
    // Step 3: Identify pivot columns
    std::vector<bool> isPivot(mat.ncol(), false);
    for (size_t i = 0; i < valid_row; i++) {
        for (size_t j = 0; j < mat.ncol(); j++) {
            if (augmented_mat(i, j) != 0) {
                isPivot[j] = true;
                break;
            }
        }
    }

    // Step 4: Solve the system for each free variable

    std::vector<std::vector<double>> null_space_vectors;
    for (size_t i = 0; i < isPivot.size(); i++) {
        if (!isPivot[i]) {
            std::vector<double> special_solution(mat.ncol(), 0);
            special_solution[i] = 1;
            for (size_t r = 0; r < valid_row; r++) {
                if (augmented_mat(r, i) != 0) {
                    special_solution[r] = -augmented_mat(r, i);
                }
            }
            null_space_vectors.push_back(special_solution);
        }
    }

    // Return the null space
    return null_space_vectors;
}

Matrix householder(std::vector<double>& x, size_t n, size_t e)
{
    // x = x - ||x||e
    x[e] = x[e] - sqrt(dot_product(x, x));

    // x/||x||
    x = normalize(x);

    Matrix H = Matrix::Identity(n, n);
    
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; j++)
        {
            H(i, j) = H(i, j) - x[i] * x[j] * 2;
        }
    }
    // h = householder matrix = I - 2 * vvt
    

    return H;
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