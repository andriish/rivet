#ifndef MATH_MATRIXN
#define MATH_MATRIXN

#include "Rivet/Math/MathHeader.hh"
#include "Rivet/Math/Vectors.hh"

//class Matrix3;
template <size_t N>
class Matrix;
template <size_t N>
Matrix<N> transpose(const Matrix<N>& m);
template <size_t N>
Matrix<N> inverse(const Matrix<N>& m);
template <size_t N>
double det(const Matrix<N>& m);
template <size_t N>
double trace(const Matrix<N>& m);
template <size_t N>
Matrix<N> multiply(const Matrix<N>& a, const Matrix<N>& b);
template <size_t N>
Matrix<N> divide(const Matrix<N>&, const double);
template <size_t N>
Matrix<N> operator*(const Matrix<N>& a, const Matrix<N>& b);

template <size_t N>
string toString(const Matrix<N>& m);
template <size_t N>
inline ostream& operator<<(std::ostream& out, const Matrix<N>& m);

template <typename Real>
bool isZero(const Real a, const Real tolerance);

typedef Matrix<4> Matrix4;


///////////////////////////////////


template <size_t N>
class Matrix {
  //friend class Matrix3;
  template <size_t M>
  friend Matrix<M> add(const Matrix<M>&, const Matrix<M>&);
  template <size_t M>
  friend Matrix<M> multiply(const double, const Matrix<M>&);
  template <size_t M>
  friend Matrix<M> multiply(const Matrix<M>&, const Matrix<M>&);
  template <size_t M>
  friend Vector<M> multiply(const Matrix<M>&, const Vector<M>&);
  template <size_t M>
  friend Matrix<M> divide(const Matrix<M>&, const double);
  template <size_t M>
  friend Matrix<M> transpose(const Matrix<M>&);
  template <size_t M>
  friend Matrix<M> inverse(const Matrix<M>&);
  template <size_t M>
  friend double det(const Matrix<M>&);
  template <size_t M>
  friend double trace(const Matrix<M>&);

public:
  static Matrix<N> Zero() {
    Matrix<N> rtn;
    return rtn;
  }

  static Matrix<N> Diag(Vector<N> diag) {
    Matrix<N> rtn;
    for (size_t i = 0; i < N; ++i) {
      rtn.set(i, i, diag[i]);
    }
    return rtn;    
  }
  
  static Matrix<N> Identity() {    
    Matrix<N> rtn;
    for (size_t i = 0; i < N; ++i) {
      rtn.set(i, i, 1);
    }
    return rtn;    
  }

public:
  Matrix() {
    _matrix.loadZero();
  }

  Matrix(const Matrix<N>& other) {
    _matrix = other._matrix;
  }
  
  Matrix& set(const size_t i, const size_t j, const double value) {
    if (i < N && j < N) {
      _matrix(i, j) = value;
    } else {
      throw std::runtime_error("Attempted set access outside matrix bounds.");
    }
    return *this;
  }
  
  double get(const size_t i, const size_t j) const {
    if (i < N && j < N) {
      return _matrix(i, j);
    } else {
      throw std::runtime_error("Attempted get access outside matrix bounds.");
    }
  }

  Vector<N> getRow(const size_t row) const {
    Vector<N> rtn;
    for (size_t i = 0; i < N; ++i) {
      rtn.set(i, _matrix(row,i));
    }
    return rtn;
  }

  Matrix<N>& setRow(const size_t row, const Vector<N>& r) {
    for (size_t i = 0; i < N; ++i) {
      _matrix(row,i) = r.get(i);
    }
    return *this;
  }

  Vector<N> getColumn(const size_t col) const {
    const Eigen::Vector<double,N> eVec = _matrix.column(col);
    Vector<N> rtn;
    for (size_t i = 0; i < N; ++i) {
      rtn.set(i, _matrix(i,col));
    }
    return rtn;
  }

  Matrix<N>& setColumn(const size_t col, const Vector<N>& c) {
    for (size_t i = 0; i < N; ++i) {
      _matrix(i,col) = c.get(i);
    }
    return *this;
  }


  Matrix<N> transpose() const {
    return ::transpose(*this);
  }
  
  Matrix<N> inverse() const {
    return ::inverse(*this);
  }

  double det() const  {
    return ::det(*this);
  }

  double trace() const {
    return ::trace(*this);
  }

  Matrix<N> operator-() const {
    Matrix<N> rtn;
    rtn._matrix = -_matrix;
    return rtn;
  }

  const size_t size() const {
    return N;
  }

  const bool isZero() const {
    for (size_t i=0; i < N; ++i) {
      for (size_t j=0; j < N; ++j) {
        if (! ::isZero(_matrix(i,j)) ) return false;
      }
    }
    return true;
  }

  const bool isEqual(Matrix<N> other) const {
    for (size_t i=0; i < N; ++i) {
      for (size_t j=i; j < N; ++j) {
        if (! ::isZero(_matrix(i,j) - other._matrix(i,j)) ) return false;
      }
    }
    return true;
  }

  const bool isSymm() const {
    return isEqual(this->transpose());
  }

  const bool isDiag() const {
    for (size_t i=0; i < N; ++i) {
      for (size_t j=0; j < N; ++j) {
        if (i == j) continue;
        if (! ::isZero(_matrix(i,j)) ) return false;
      }
    }
    return true;
  }

  bool operator==(const Matrix<N>& a) const {
    return _matrix == a._matrix;
  }

  bool operator!=(const Matrix<N>& a) const {
    return _matrix != a._matrix;
  }

  bool operator<(const Matrix<N>& a) const {
    return _matrix < a._matrix;
  }

  bool operator<=(const Matrix<N>& a) const {
    return _matrix <= a._matrix;
  }

  bool operator>(const Matrix<N>& a) const {
    return _matrix > a._matrix;
  }

  bool operator>=(const Matrix<N>& a) const {
    return _matrix >= a._matrix;
  }

  Matrix<N>& operator*=(const Matrix<N>& m) {
    _matrix = _matrix * m._matrix;
    return *this;
  }

  Matrix<N>& operator*=(const double a) {
    _matrix *= a;
    return *this;
  }

  Matrix<N>& operator/=(const double a) {
    _matrix /= a;
    return *this;
  }

  Matrix<N>& operator+=(const Matrix<N>& m) {
    _matrix += m._matrix;
    return *this;
  }

  Matrix<N>& operator-=(const Matrix<N>& m) {
    _matrix -= m._matrix;
    return *this;
  }

protected:
  typedef Eigen::Matrix<double,N> EMatrix;
  EMatrix _matrix;
};


/////////////////////////////////


template <size_t N>
inline Matrix<N> add(const Matrix<N>& a, const Matrix<N>& b) {
  Matrix<N> result;
  result._matrix = a._matrix + b._matrix;
  return result;
}

template <size_t N>
inline Matrix<N> subtract(const Matrix<N>& a, const Matrix<N>& b) {
  return add(a, -b);
}

template <size_t N>
inline Matrix<N> operator+(const Matrix<N> a, const Matrix<N>& b) {
  return add(a, b);
}

template <size_t N>
inline Matrix<N> operator-(const Matrix<N> a, const Matrix<N>& b) {
  return subtract(a, b);
}

template <size_t N>
inline Matrix<N> multiply(const double a, const Matrix<N>& m) {
  Matrix<N> rtn;
  rtn._matrix = a * m._matrix;
  return rtn;
}

template <size_t N>
inline Matrix<N> multiply(const Matrix<N>& m, const double a) {
  return multiply(a, m);
}

template <size_t N>
inline Matrix<N> divide(const Matrix<N>& m, const double a) {
  return multiply(1/a, m);
}

template <size_t N>
inline Matrix<N> operator*(const double a, const Matrix<N>& m) {
  return multiply(a, m);
}

template <size_t N>
inline Matrix<N> operator*(const Matrix<N>& m, const double a) {
  return multiply(a, m);
}

template <size_t N>
inline Matrix<N> multiply(const Matrix<N>& a, const Matrix<N>& b) {
  Matrix<N> tmp;
  tmp._matrix = a._matrix * b._matrix;
  return tmp;
}

template <size_t N>
inline Matrix<N> operator*(const Matrix<N>& a, const Matrix<N>& b) {
  return multiply(a, b);
}


template <size_t N>
inline Vector<N> multiply(const Matrix<N>& a, const Vector<N>& b) {
  Vector<N> tmp;
  tmp._vec = a._matrix * b._vec;
  return tmp;
}

template <size_t N>
inline Vector<N> operator*(const Matrix<N>& a, const Vector<N>& b) {
  return multiply(a, b);
}

template <size_t N>
inline Matrix<N> transpose(const Matrix<N>& m) {
  Matrix<N> tmp;
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      tmp.set(i, j, m.get(j, i));
    }
  }
  return tmp;
}

template <size_t N>
inline Matrix<N> inverse(const Matrix<N>& m) {
  Matrix<N> tmp;
  tmp._matrix = m._matrix.inverse();
  return tmp;
}

template <size_t N>
inline double det(const Matrix<N>& m) {
  return m._matrix.determinant();
}

template <size_t N>
inline double trace(const Matrix<N>& m) {
  double tr = 0;
  for (size_t i = 0; i < N; ++i) {
    tr += m._matrix(i,i);
  }
  return tr;
  //m._matrix.trace();
}


/////////////////////////////////


template <size_t N>
inline string toString(const Matrix<N>& m) {
  ostringstream ss;
  ss << "[ ";
  for (size_t i = 0; i < m.size(); ++i) {
    ss << "( ";
    for (size_t j = 0; j < m.size(); ++j) {
      const double e = m.get(i, j);
      ss << (isZero(e) ? 0.0 : e) << " ";
    }
    ss << ") ";
  }
  ss << "]";
  return ss.str();
}


template <size_t N>
inline ostream& operator<<(std::ostream& out, const Matrix<N>& m) {
  out << toString(m);
  return out;
}


#endif
