// Tell swig the name of the module we're creating
%module pymass

// Pull in the headers from Python itself and from our library
%{
#define SWIG_FILE_WITH_INIT
#define SWIG_PYTHON_STRICT_BYTE_CHAR
#include <Python.h>
#include "pymass_export.h"
#include "MZXML.h"
#include "utils.h"
%}

%include <typemaps.i>
%include <std_string.i>
%include "std_vector.i"


namespace std {
   %template(IntVector) vector<int>;
   %template(DoubleVector) vector<double>;
   %template(StringVector) vector<string>;
   %template(ConstCharVector) vector<const char*>;
}


// Eigen matrices into Numpy arrays.
%include <eigen.i>
%eigen_typemaps(Eigen::VectorXd)
%eigen_typemaps(Eigen::MatrixXd)
%eigen_typemaps(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>)



%include <windows.i>
// Tell swig to build bindings for everything in our library
%include "pymass_export.h"
%include "MZXML.h"
%include "utils.h"

