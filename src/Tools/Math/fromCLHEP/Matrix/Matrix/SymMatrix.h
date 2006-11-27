// -*- C++ -*-
// CLASSDOC OFF
// $Id: SymMatrix.h,v 1.4 2005/04/27 19:31:54 garren Exp $
// ---------------------------------------------------------------------------
// CLASSDOC ON
// 
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
// 
// This is the definition of the HepSymMatrix class.
// 
// Copyright Cornell University 1993, 1996, All Rights Reserved.
// 
// This software written by Nobu Katayama and Mike Smyth, Cornell University.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 1. Redistributions of source code must retain the above copyright
//    notice and author attribution, this list of conditions and the
//    following disclaimer. 
// 2. Redistributions in binary form must reproduce the above copyright
//    notice and author attribution, this list of conditions and the
//    following disclaimer in the documentation and/or other materials
//    provided with the distribution.
// 3. Neither the name of the University nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// 
// Creation of derivative forms of this software for commercial
// utilization may be subject to restriction; written permission may be
// obtained from Cornell University.
// 
// CORNELL MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.  By way
// of example, but not limitation, CORNELL MAKES NO REPRESENTATIONS OR
// WARRANTIES OF MERCANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT
// THE USE OF THIS SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS,
// COPYRIGHTS, TRADEMARKS, OR OTHER RIGHTS.  Cornell University shall not be
// held liable for any liability with respect to any claim by the user or any
// other party arising from use of the program.
//
// .SS Usage
//
//   This is very much like the Matrix, except of course it is restricted to
//   Symmetric Matrix.  All the operations for Matrix can also be done here
//   (except for the +=,-=,*= that don't yield a symmetric matrix.  e.g.
//    +=(const Matrix &) is not defined)
   
//   The matrix is stored as a lower triangular matrix.
//   In addition to the (row, col) method of finding element, fast(row, col)
//   returns an element with checking to see if row and col need to be 
//   interchanged so that row >= col.

//   New operations are:
//
// .ft B
//  sym = s.similarity(m);
//
//  This returns m*s*m.T(). This is a similarity
//  transform when m is orthogonal, but nothing
//  restricts m to be orthogonal.  It is just
//  a more efficient way to calculate m*s*m.T,
//  and it realizes that this should be a 
//  HepSymMatrix (the explicit operation m*s*m.T
//  will return a Matrix, not realizing that 
//  it is symmetric).
//
// .ft B
//  sym =  similarityT(m);
//
// This returns m.T()*s*m.
//
// .ft B
// s << m;
//
// This takes the matrix m, and treats it
// as symmetric matrix that is copied to s.
// This is useful for operations that yield
// symmetric matrix, but which the computer
// is too dumb to realize.
//
// .ft B
// s = vT_times_v(const HepVector v);
//
//  calculates v.T()*v.
//
// ./"This code has been written by Mike Smyth, and the algorithms used are
// ./"described in the thesis "A Tracking Library for a Silicon Vertex Detector"
// ./"(Mike Smyth, Cornell University, June 1993).
// ./"Copyright (C) Cornell University 1993. Permission is granted to copy and 
// ./"distribute this code, provided this copyright is not changed or deleted.
// ./"You may modify your copy, providing that you cause the modified file to
// ./"carry prominent notices stating that you changed the files, and the date
// ./"of any change. This code may not be sold, nor may it be contained in
// ./"programs that are to be sold without the written permission of the author.
// ./"You may, however, charge a fee for the physical act of transferring a
// ./"copy of this code. The code is offered "as is" without warranty of any 
// ./"kind, either expressed or implied.  By copying, distributing, or 
// ./"modifying this code you indicate your acceptance of this license to
// ./"do so, and all its terms and conditions.
// ./"This is file contains C++ stuff for doing things with Matrixes.
// ./"To turn on bound checking, define MATRIX_BOUND_CHECK before including
// ./"this file.
//

#ifndef _SYMMatrix_H_
#define _SYMMatrix_H_

#ifdef GNUPRAGMA
#pragma interface
#endif

#include <vector>

#include "CLHEP/Matrix/defs.h"
#include "CLHEP/Matrix/GenMatrix.h"

namespace CLHEP {

class HepRandom;

class HepMatrix;
class HepDiagMatrix;
class HepVector;

/**
 * @author
 * @ingroup matrix
 */
class HepSymMatrix : public HepGenMatrix {
public:
   inline HepSymMatrix();
   // Default constructor. Gives 0x0 symmetric matrix.
   // Another SymMatrix can be assigned to it.

   explicit HepSymMatrix(int p);
   HepSymMatrix(int p, int);
   // Constructor. Gives p x p symmetric matrix.
   // With a second argument, the matrix is initialized. 0 means a zero
   // matrix, 1 means the identity matrix.

   HepSymMatrix(int p, HepRandom &r);

   HepSymMatrix(const HepSymMatrix &m1);
   // Copy constructor.

   HepSymMatrix(const HepDiagMatrix &m1);
   // Constructor from DiagMatrix

   virtual ~HepSymMatrix();
   // Destructor.

   inline int num_row() const;
   inline int num_col() const;
   // Returns number of rows/columns.

   const double & operator()(int row, int col) const; 
   double & operator()(int row, int col);
   // Read and write a SymMatrix element.
   // ** Note that indexing starts from (1,1). **

   const double & fast(int row, int col) const;
   double & fast(int row, int col);
   // fast element access.
   // Must be row>=col;
   // ** Note that indexing starts from (1,1). **

   void assign(const HepMatrix &m2);
   // Assigns m2 to s, assuming m2 is a symmetric matrix.

   void assign(const HepSymMatrix &m2);
   // Another form of assignment. For consistency.

   HepSymMatrix & operator*=(double t);
   // Multiply a SymMatrix by a floating number.

   HepSymMatrix & operator/=(double t); 
   // Divide a SymMatrix by a floating number.

   HepSymMatrix & operator+=( const HepSymMatrix &m2);
   HepSymMatrix & operator+=( const HepDiagMatrix &m2);
   HepSymMatrix & operator-=( const HepSymMatrix &m2);
   HepSymMatrix & operator-=( const HepDiagMatrix &m2);
   // Add or subtract a SymMatrix.

   HepSymMatrix & operator=( const HepSymMatrix &m2);
   HepSymMatrix & operator=( const HepDiagMatrix &m2);
   // Assignment operators. Notice that there is no SymMatrix = Matrix.

   HepSymMatrix operator- () const;
   // unary minus, ie. flip the sign of each element.

   HepSymMatrix T() const;
   // Returns the transpose of a SymMatrix (which is itself).

   HepSymMatrix apply(double (*f)(double, int, int)) const;
   // Apply a function to all elements of the matrix.

   HepSymMatrix similarity(const HepMatrix &m1) const;
   HepSymMatrix similarity(const HepSymMatrix &m1) const;
   // Returns m1*s*m1.T().

   HepSymMatrix similarityT(const HepMatrix &m1) const;
   // temporary. test of new similarity.
   // Returns m1.T()*s*m1.

   double similarity(const HepVector &v) const;
   // Returns v.T()*s*v (This is a scaler).

   HepSymMatrix sub(int min_row, int max_row) const;
   // Returns a sub matrix of a SymMatrix.
   void sub(int row, const HepSymMatrix &m1);
   // Sub matrix of this SymMatrix is replaced with m1.
   HepSymMatrix sub(int min_row, int max_row);
   // SGI CC bug. I have to have both with/without const. I should not need
   // one without const.

   inline HepSymMatrix inverse(int &ifail) const;
   // Invert a Matrix. The matrix is not changed
   // Returns 0 when successful, otherwise non-zero.

   void invert(int &ifail);
   // Invert a Matrix.
   // N.B. the contents of the matrix are replaced by the inverse.
   // Returns ierr = 0 when successful, otherwise non-zero. 
   // This method has less overhead then inverse().

   double determinant() const;
   // calculate the determinant of the matrix.

   double trace() const;
   // calculate the trace of the matrix (sum of diagonal elements).

   class HepSymMatrix_row {
   public:
      inline HepSymMatrix_row(HepSymMatrix&,int);
      inline double & operator[](int);
   private:
      HepSymMatrix& _a;
      int _r;
   };
   class HepSymMatrix_row_const {
   public:
      inline HepSymMatrix_row_const(const HepSymMatrix&,int);
      inline const double & operator[](int) const;
   private:
      const HepSymMatrix& _a;
      int _r;
   };
   // helper class to implement m[i][j]

   inline HepSymMatrix_row operator[] (int);
   inline HepSymMatrix_row_const operator[] (int) const;
   // Read or write a matrix element.
   // While it may not look like it, you simply do m[i][j] to get an
   // element. 
   // ** Note that the indexing starts from [0][0]. **

   // Special-case inversions for 5x5 and 6x6 symmetric positive definite:
   // These set ifail=0 and invert if the matrix was positive definite;
   // otherwise ifail=1 and the matrix is left unaltered.
   void invertCholesky5 (int &ifail);  
   void invertCholesky6 (int &ifail);

   // Inversions for 5x5 and 6x6 forcing use of specific methods:  The
   // behavior (though not the speed) will be identical to invert(ifail).
   void invertHaywood4 (int & ifail);  
   void invertHaywood5 (int &ifail);  
   void invertHaywood6 (int &ifail);
   void invertBunchKaufman (int &ifail);  

protected:
   inline int num_size() const;
  
private:
   friend class HepSymMatrix_row;
   friend class HepSymMatrix_row_const;
   friend class HepMatrix;
   friend class HepDiagMatrix;

   friend void tridiagonal(HepSymMatrix *a,HepMatrix *hsm);
   friend double condition(const HepSymMatrix &m);
   friend void diag_step(HepSymMatrix *t,int begin,int end);
   friend void diag_step(HepSymMatrix *t,HepMatrix *u,int begin,int end);
   friend HepMatrix diagonalize(HepSymMatrix *s);
   friend HepVector house(const HepSymMatrix &a,int row,int col);
   friend void house_with_update2(HepSymMatrix *a,HepMatrix *v,int row,int col);

   friend HepSymMatrix operator+(const HepSymMatrix &m1, 
				  const HepSymMatrix &m2);
   friend HepSymMatrix operator-(const HepSymMatrix &m1, 
				  const HepSymMatrix &m2);
   friend HepMatrix operator*(const HepSymMatrix &m1, const HepSymMatrix &m2);
   friend HepMatrix operator*(const HepSymMatrix &m1, const HepMatrix &m2);
   friend HepMatrix operator*(const HepMatrix &m1, const HepSymMatrix &m2);
   friend HepVector operator*(const HepSymMatrix &m1, const HepVector &m2);
   // Multiply a Matrix by a Matrix or Vector.
   
   friend HepSymMatrix vT_times_v(const HepVector &v);
   // Returns v * v.T();

#ifdef DISABLE_ALLOC
   std::vector<double > m;
#else
   std::vector<double,Alloc<double,25> > m;
#endif
   int nrow;
   int size;				     // total number of elements

   static double posDefFraction5x5;
   static double adjustment5x5;
   static const  double CHOLESKY_THRESHOLD_5x5;
   static const  double CHOLESKY_CREEP_5x5;

   static double posDefFraction6x6;
   static double adjustment6x6;
   static const double CHOLESKY_THRESHOLD_6x6;
   static const double CHOLESKY_CREEP_6x6;

   void invert4  (int & ifail);
   void invert5  (int & ifail);
   void invert6  (int & ifail);

};

//
// Operations other than member functions for Matrix, SymMatrix, DiagMatrix
// and Vectors implemented in Matrix.cc and Matrix.icc (inline).
//

std::ostream& operator<<(std::ostream &s, const HepSymMatrix &q);
// Write out Matrix, SymMatrix, DiagMatrix and Vector into ostream.

HepMatrix operator*(const HepMatrix &m1, const HepSymMatrix &m2);
HepMatrix operator*(const HepSymMatrix &m1, const HepMatrix &m2);
HepMatrix operator*(const HepSymMatrix &m1, const HepSymMatrix &m2);
HepSymMatrix operator*(double t, const HepSymMatrix &s1);
HepSymMatrix operator*(const HepSymMatrix &s1, double t);
// Multiplication operators.
// Note that m *= m1 is always faster than m = m * m1

HepSymMatrix operator/(const HepSymMatrix &m1, double t);
// s = s1 / t. (s /= t is faster if you can use it.)

HepMatrix operator+(const HepMatrix &m1, const HepSymMatrix &s2);
HepMatrix operator+(const HepSymMatrix &s1, const HepMatrix &m2);
HepSymMatrix operator+(const HepSymMatrix &s1, const HepSymMatrix &s2);
// Addition operators

HepMatrix operator-(const HepMatrix &m1, const HepSymMatrix &s2);
HepMatrix operator-(const HepSymMatrix &m1, const HepMatrix &m2);
HepSymMatrix operator-(const HepSymMatrix &s1, const HepSymMatrix &s2);
// subtraction operators

HepSymMatrix dsum(const HepSymMatrix &s1, const HepSymMatrix &s2);
// Direct sum of two symmetric matrices;

double condition(const HepSymMatrix &m);
// Find the conditon number of a symmetric matrix.

void diag_step(HepSymMatrix *t, int begin, int end);
void diag_step(HepSymMatrix *t, HepMatrix *u, int begin, int end);
// Implicit symmetric QR step with Wilkinson Shift

HepMatrix diagonalize(HepSymMatrix *s);
// Diagonalize a symmetric matrix.
// It returns the matrix U so that s_old = U * s_diag * U.T()

HepVector house(const HepSymMatrix &a, int row=1, int col=1);
void house_with_update2(HepSymMatrix *a, HepMatrix *v, int row=1, int col=1);
// Finds and does Householder reflection on matrix.

void tridiagonal(HepSymMatrix *a, HepMatrix *hsm);
HepMatrix tridiagonal(HepSymMatrix *a);
// Does a Householder tridiagonalization of a symmetric matrix.

}  // namespace CLHEP

#ifdef ENABLE_BACKWARDS_COMPATIBILITY
//  backwards compatibility will be enabled ONLY in CLHEP 1.9
using namespace CLHEP;
#endif

#ifndef HEP_DEBUG_INLINE
#include "CLHEP/Matrix/SymMatrix.icc"
#endif

#endif /*!_SYMMatrix_H*/
