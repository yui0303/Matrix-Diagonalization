#include "_matrix.hpp"
#include "_eigen.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sstream>

namespace py = pybind11;

PYBIND11_MODULE(_matrix, m) {
    py::class_<Matrix>(m, "Matrix")
        .def(py::init<size_t, size_t>())
        .def(py::init<size_t, size_t, std::vector<double> const &>())
        .def(py::init<Matrix const&>())
        .def("__getitem__", [](Matrix &mat, std::pair<size_t, size_t> index){
            return mat(index.first, index.second);
        })
        .def("__setitem__", [](Matrix &mat, std::pair<size_t, size_t> index, double val){
            mat(index.first, index.second) = val;
        })
        .def("__eq__", &Matrix::operator==)
        .def("__repr__", [](const Matrix &mat){
            std::ostringstream os;
            os << mat;
            return os.str();
        }) 
        .def("__mul__", [](Matrix &mat1, Matrix &mat2){
            return mat1 * mat2;
        })
        .def("transpose", &Matrix::transpose)
        .def_property_readonly("nrow", &Matrix::nrow)
        .def_property_readonly("ncol", &Matrix::ncol)
        ;
    m.def("multiply_tile", &multiply_tile);
    m.def("multiply_naive", &multiply_naive);
    m.def("QR_decomposition", &QR_decomposition_GS);
    // m.def("multiply_mkl", &multiply_mkl);
}