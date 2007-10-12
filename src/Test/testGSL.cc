#include <iostream>
#include <sstream>
#include <cmath>
#include <string>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_eigen.h"

using namespace std;


void gsl_symm_matrix_set(gsl_matrix* A, size_t i, size_t j, double val) {
  gsl_matrix_set(A, i, j, val);
  if (i != j) {
    gsl_matrix_set(A, j, i, val);
  }
}

const string toString(gsl_vector* v) {
  ostringstream ss;
  ss << "( ";
  for (size_t i = 0; i < v->size; ++i) {
    const double e = gsl_vector_get(v, i);
    ss << e << " ";
  }
  ss << ")";
  return ss.str();
}

const string toString(gsl_matrix* m) {
  ostringstream ss;
  ss << "[ ";
  for (size_t i = 0; i < m->size1; ++i) {
    ss << "( ";
    for (size_t j = 0; j < m->size2; ++j) {
      const double e = gsl_matrix_get(m, i, j);
      ss << e << " ";
    }
    ss << ") ";
  }
  ss << "]";
  return ss.str();
}


int main() {
  gsl_matrix* A = gsl_matrix_alloc(3, 3);
  gsl_matrix_set_zero(A);
  gsl_symm_matrix_set(A, 0, 0,  7/4.0);
  gsl_symm_matrix_set(A, 1, 1, 13/4.0);
  gsl_symm_matrix_set(A, 2, 2,  9);
  gsl_symm_matrix_set(A, 0, 1, 3*sqrt(3)/4.0);

  cout << toString(A) << endl;

  gsl_matrix* vecs = gsl_matrix_alloc(3, 3);
  gsl_vector* vals = gsl_vector_alloc(3);
  gsl_eigen_symmv_workspace* workspace = gsl_eigen_symmv_alloc (3);
  const int status = gsl_eigen_symmv(A, vals, vecs, workspace);

  cout << "evals = " << toString(vals) << endl; 
  cout << "evecs = " << toString(vecs) << endl;

  gsl_eigen_symmv_free(workspace);
  gsl_matrix_free(A);
  gsl_matrix_free(vecs);
  gsl_vector_free(vals);

  return EXIT_SUCCESS;
}
