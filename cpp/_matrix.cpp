#include "_matrix.hpp"

Matrix::Matrix(size_t nrow, size_t ncol) : m_nrow(nrow), m_ncol(ncol)
{
    size_t nelement = nrow * ncol;
    m_buffer = new double[nelement];
    for (size_t i=0; i<nelement; ++i)
    {
        m_buffer[i] = 0;
    }
}

Matrix::Matrix(size_t nrow, size_t ncol, std::vector<double> const & vec) : m_nrow(nrow), m_ncol(ncol)
{
    if (vec.size() != nrow * ncol)
    {
        throw std::out_of_range("Matrix::Matrix(): vector size differs from matrix size");
    }

    m_buffer = new double[nrow * ncol];

    for (size_t i=0; i<nrow; ++i)
    {
        for (size_t j=0; j<ncol; ++j)
        {
            m_buffer[i*ncol + j] = vec[i*ncol + j];
        }
    }

}

Matrix::Matrix(Matrix const &mat) : m_nrow(mat.m_nrow), m_ncol(mat.m_ncol)
{
    m_buffer = new double[m_nrow * m_ncol];
    std::copy(mat.m_buffer, mat.m_buffer + m_nrow * m_ncol, m_buffer);
}

Matrix & Matrix::operator=(Matrix const & mat)
{
    m_nrow = mat.nrow();
    m_ncol = mat.ncol();
    m_buffer = new double[m_nrow * m_ncol];

    for (size_t i=0; i<m_nrow; ++i)
    {
        for (size_t j=0; j<m_ncol; ++j)
        {
            m_buffer[i*m_ncol + j] = mat(i,j);
        }
    }
    return *this;
}

bool Matrix::operator==(Matrix const & mat) const
{
    if (mat.m_ncol != m_ncol || mat.m_nrow != m_nrow) return false;
    for (size_t i = 0; i < mat.m_nrow; ++i)
    {
        for (size_t j = 0; j < mat.m_ncol; ++j)
        {
            if(mat(i, j) != m_buffer[i*m_nrow+j])return false;
        }
    }
    return true;
}

std::ostream & operator<<(std::ostream & os, Matrix const &mat)
{
    for (size_t i = 0; i < mat.m_nrow; ++i)
    {
        for (size_t j = 0; j < mat.m_ncol; ++j)
        {
            os << std::setw(10) << std::fixed << std::setprecision(8) << std::right << mat(i, j) << " ";
        }
        os << std::endl;
    }
    return os;
}


double   Matrix::operator() (size_t row, size_t col) const
{
    //add the boundary check
    if (row >= m_nrow || col >= m_ncol)
    {
        throw std::out_of_range("Matrix::operator(): index out of range");
    }

    return m_buffer[row*m_ncol + col];
}

double & Matrix::operator() (size_t row, size_t col)
{
    //add the boundary check
    if (row >= m_nrow || col >= m_ncol)
    {
        throw std::out_of_range("Matrix::operator(): index out of range");
    }

    return m_buffer[row*m_ncol + col];
}

// return a row vector
std::vector<double> Matrix::operator() (size_t row) const
{
    if (row >= m_nrow)
    {
        throw std::out_of_range("Matrix::operator(): index out of range");
    }

    std::vector<double> vec(m_ncol);
    for (size_t i=0; i<m_ncol; ++i)
    {
        vec[i] = m_buffer[row*m_ncol + i];
    }
    return vec;
}

Matrix & Matrix::operator*(Matrix const & mat)
{
    if (m_ncol != mat.m_nrow)
    {
        throw std::out_of_range("Matrix::operator*: invalid matrix dimensions for multiplication");
    }

    Matrix *result = new Matrix(m_nrow, mat.m_ncol);

    for (size_t i=0; i<m_nrow; ++i)
    {
        for (size_t j=0; j<mat.m_ncol; ++j)
        {
            double v = 0;
            for (size_t k=0; k<m_ncol; ++k)
            {
                v += (*this)(i,k) * mat(k,j);
            }
            (*result)(i, j) = v;
        }
    }

    return *result;
}

Matrix & Matrix::operator*(std::vector<double> const & vec)
{
    if (m_ncol != vec.size())
    {
        throw std::out_of_range("Matrix::operator*: matrix column differs from vector size");
    }

    Matrix *result = new Matrix(m_nrow, 1);

    for (size_t i=0; i<m_nrow; ++i)
    {
        double v = 0;
        for (size_t j=0; j<m_ncol; ++j)
        {
            v += (*this)(i,j) * vec[j];
        }
        (*result)(i, 0) = v;
    }

    return *result;
}

Matrix Matrix::transpose() const
{
    Matrix result(m_ncol, m_nrow);

    for (size_t i=0; i<m_nrow; ++i)
    {
        for (size_t j=0; j<m_ncol; ++j)
        {
            result(j, i) = (*this)(i, j);
        }
    }

    return result;
}

double* Matrix::data() const { return m_buffer; }

size_t Matrix::nrow() const { return m_nrow; }
size_t Matrix::ncol() const { return m_ncol; }

Matrix::~Matrix()
{
    delete[] m_buffer;
}

//operator overloading for vector
std::ostream & operator<<(std::ostream & os, std::vector<double> const & vec)
{
    for (size_t i=0; i<vec.size(); ++i)
    {
        os << vec[i] << " ";
    }
    os << std::endl;
    return os;
}
