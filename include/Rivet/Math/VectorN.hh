#ifndef RIVET_MATH_VECTORN
#define RIVET_MATH_VECTORN

#include "Rivet/Math/MathHeader.hh"

template <size_t N>
class Vector;
template <size_t N>
class Matrix;

template <size_t N>
double mod2(const Vector<N>& v);
template <size_t N>
double mod(const Vector<N>& v);
template <size_t N>
Vector<N> multiply(const Matrix<N>& a, const Vector<N>& b);

template <typename Real>
bool isZero(const Real a, const Real tolerance);


/// A minimal base class for ND vectors.
template <size_t N>
class Vector {
  template <size_t M>
  friend Vector<M> multiply(const Matrix<M>& a, const Vector<M>& b);

public:
  Vector() { _vec.loadZero(); }

  Vector(const Vector<N>& other) 
    : _vec(other._vec) { }

  const double& get(const size_t index) const {
    if (index >= N) {
      throw std::runtime_error("Tried to access an invalid vector index.");
    } else {
      return _vec(index);
    }
  }

  /// Direct access to vector elements by index.
  const double& operator[](const size_t index) const {
    return get(index);
  }

  /// Direct access to vector elements by index.
  double& operator[](const size_t index) {
    return get(index);
  }

  Vector<N>& set(const size_t index, const double value) {
    if (index >= N) {
      throw std::runtime_error("Tried to access an invalid vector index.");
    } else {
      _vec[index] = value;
    }
    return *this;
  }

  const size_t size() const {
    return N;
  }

  const bool isZero() const {
    for (size_t i=0; i < N; ++i) {
      if (! ::isZero(_vec[i]) ) return false;
    }
    return true;
  }

  Vector<N> operator-() const {
    Vector<N> rtn;
    rtn._vec = -_vec;
    return rtn;
  }

  inline bool operator==(const Vector<N>& a) {
    return _vec == a._vec;
  }

  inline bool operator!=(const Vector<N>& a) {
    return _vec != a._vec;
  }

  inline bool operator<(const Vector<N>& a) {
    return _vec < a._vec;
  }

  inline bool operator<=(const Vector<N>& a) {
    return _vec <= a._vec;
  }

  inline bool operator>(const Vector<N>& a) {
    return _vec > a._vec;
  }

  inline bool operator>=(const Vector<N>& a) {
    return _vec >= a._vec;
  }


protected:
  double& get(const size_t index) {
    if (index >= N) {
      throw std::runtime_error("Tried to access an invalid vector index.");
    } else {
      return _vec(index);
    }
  }

  /// Vector
  Eigen::Vector<double,N> _vec;
};



/// Calculate the modulus-squared of a vector.
/// \f$ \sum_{i=1}^N x_i^2 \f$.
template <size_t N>
inline double mod2(const Vector<N>& v) {
  double mod2 = 0.0;
  for (size_t i = 0; i < v.size(); ++i) {
    const double element = v.get(i);
    mod2 += element*element;
  }
  return mod2;
}

/// Calculate the modulus of a vector.
/// \f$ \sqrt{\sum_{i=1}^N x_i^2} \f$.
template <size_t N>
inline double mod(const Vector<N>& v) {
  const double norm = mod2(v);
  assert(norm >= 0);
  return sqrt(norm);
}


/////////////////////////////////////////////////

template <size_t N>
inline const string toString(const Vector<N>& v) {
  ostringstream out;
  out << "(";
  for (size_t i = 0; i < v.size(); ++i) {
    out << (fabs(v[i]) < 1E-30 ? 0.0 : v[i]);
    if (i < v.size()-1) out << ", ";
  }
  out << ")";
  return out.str();
}

template <size_t N>
inline std::ostream& operator<<(std::ostream& out, const Vector<N>& v) {
  out << toString(v);
  return out;
}

#endif
