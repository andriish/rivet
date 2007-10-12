#include <stdexcept>
#include <string>
#include <ostream>
#include <sstream>
#include <iostream>
#include <cmath>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_eigen.h"


using std::string;
using std::ostringstream;
using std::cout;
using std::endl;


//template <int N>
//class Matrix<N>;




template <int N>
class Matrix {
public:
  Matrix() : _size(N) {
    for (size_t i = 0; i < _size; ++i) {
      for (size_t j = 0; j < _size; ++j) {
        set(i, j, 0);
      }
    }
  }
  
  Matrix& set(size_t i, size_t j, double val) {
    if (i < N && j < N) {
      _elements[i][j] = val;
    } else {
      throw std::runtime_error("Attempted set access outside matrix bounds.");
    }
    return *this;
  }

  double get(size_t i, size_t j) const {
    if (i < N && j < N) {
      return _elements[i][j];
    } else {
      throw std::runtime_error("Attempted get access outside matrix bounds.");
    }
  }

private:  
  double _elements[N][N];
  const size_t _size;
};


template <int N>
string toString(const Matrix<N>& m) {
  ostringstream ss;
  ss << "[ ";
  for (size_t i = 0; i < N; ++i) {
    ss << "( ";
    for (size_t j = 0; j < N; ++j) {
      const double e = m.get(i, j);
      ss << e << " ";
    }
    ss << ") ";
  }
  ss << "]";
  return ss.str();
}


  

template <int N>
class Eigensystem {
public:
  Eigensystem& solve(const Matrix<N>& m) {
    gsl_matrix* A = gsl_matrix_alloc(N, N);
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        gsl_matrix_set(A, i, j, m.get(i, j));
      }
    }

    gsl_matrix* vecs = gsl_matrix_alloc(3, 3);
    gsl_vector* vals = gsl_vector_alloc(3);
    gsl_eigen_symmv_workspace* workspace = gsl_eigen_symmv_alloc (3);
    gsl_eigen_symmv(A, vals, vecs, workspace);

    for (size_t i = 0; i < N; ++i) {
      _evals[i] = gsl_vector_get(vals, i);
      for (size_t j = 0; j < N; ++j) {
        _evecs[i][j] = gsl_matrix_get(vecs, i, j);
      }
    }

    gsl_eigen_symmv_free(workspace);
    gsl_matrix_free(A);
    gsl_matrix_free(vecs);
    gsl_vector_free(vals);

    return *this;
  }

  void show() {
    cout << "Eigenvalues = ";
    for (size_t i = 0; i < N; ++i) {
      cout << _evals[i] << " ";
    }
    cout << endl;

    cout << "Eigenvectors = ";
    for (size_t i = 0; i < N; ++i) {
      cout << "( ";
      for (size_t j = 0; j < N; ++j) {
        cout << _evecs[i][j] << " ";
      }
      cout << ") ";
    }
    cout << endl;
  }

private:
  double _evals[N];
  double _evecs[N][N];
};



class LorentzTransform : public Matrix<4> {
  
};
