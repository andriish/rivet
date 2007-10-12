#include <stdexcept>
#include <string>
#include <ostream>
#include <sstream>
#include <iostream>
#include <cmath>

using std::string;
using std::ostringstream;
using std::cout;
using std::endl;

// Forward declarations.
class Vector3;
class Vector4;
double mod2(const Vector3& v3);
double mod(const Vector3& v3);
double invariant(const Vector4& v4);



/////////////////////////////////////////////////////////////////////



/// A minimal base class for 3D vectors.
class Vector3 {
public:

  //  virtual ~Vector3() = 0;

public:
  /// Direct access to vector elements by index.
  const double& operator[](size_t index) const {
    if (index < 0 || index > 3) {
      throw std::runtime_error("Tried to access an invalid 3-vector index.");
    } else {
      return _elements[index];
    }
  }

  /// Direct access to vector elements by index.
  /// @todo Protect caching
  double& operator[](size_t index) {
    if (index < 0 || index > 3) {
      throw std::runtime_error("Tried to access an invalid 3-vector index.");
    } else {
      return _elements[index];
    }
  }


  double mod2() const {
    if (_dirty) {
      _mod2 = ::mod2(*this);
      _dirty = false;
    }
    return _mod2;
  }


  double mod() const { 
    if (_dirty) {
      _mod = ::mod(*this);
      _dirty = false;
    }
    return _mod;
  }

protected:
  /// Flag used to determine if the cached variables can be used
  mutable bool _dirty;

  /// Caching variables
  mutable double _mod2, _mod;

  /// Vector elements
  double _elements[4];
};



class Vector3X : public Vector3 {
public:
  Vector3X(double x, double y, double z) {
    this->x(x);
    this->y(y);
    this->z(z);
  }

  ~Vector3X() {}

public:
  double x() const { return _elements[0]; }
  double y() const { return _elements[1]; }
  double z() const { return _elements[2]; }
  Vector3X& x(const double x) { _elements[0] = x; _dirty = true; return *this; }
  Vector3X& y(const double y) { _elements[1] = y; _dirty = true; return *this; }
  Vector3X& z(const double z) { _elements[2] = z; _dirty = true; return *this; }
};



class Vector3P : public Vector3 {
public:
  Vector3P(double px, double py, double pz) {
    this->px(px);
    this->py(py);
    this->pz(pz);
  }

  ~Vector3P() {}

public:
  double px() const { return _elements[0]; }
  double py() const { return _elements[1]; }
  double pz() const { return _elements[2]; }
  Vector3P& px(const double px) { _elements[0] = px; _dirty = true; return *this; }
  Vector3P& py(const double py) { _elements[1] = py; _dirty = true; return *this; }
  Vector3P& pz(const double pz) { _elements[2] = pz; _dirty = true; return *this; }
};


///////////////////////////////////////////////////////////////////


class Vector4 {
public:

  //  virtual ~Vector4() = 0;

public:
  const double& operator[](size_t index) const {
    if (index < 0 || index > 3) {
      throw std::runtime_error("Tried to access an invalid 4-vector index.");
    } else {
      return _elements[index];
    }
  }

  double& operator[](size_t index)  {
    if (index < 0 || index > 3) {
      throw std::runtime_error("Tried to access an invalid 4-vector index.");
    } else {
      return _elements[index];
    }
  }
  
  double invariant() const {
    if (_dirty) {
      _invariant = ::invariant(*this);
      _dirty = false;
    } 
    return _invariant;
  }

  /// Get the spatial part of the 4-vector as a 3-vector.
  //virtual Vector3 vector3() const = 0;

protected:
  mutable bool _dirty;

  mutable double _invariant;

  double _elements[4];
};


class Vector4X : public Vector4 {
public:
  Vector4X(double t, double x, double y, double z) {
    this->t(t);
    this->x(x);
    this->y(y);
    this->z(z);
  }

  ~Vector4X() {}

public:
  double t() const { return _elements[0]; }
  double x() const { return _elements[1]; }
  double y() const { return _elements[2]; }
  double z() const { return _elements[3]; }
  Vector4X& t(double t) { _elements[0] = t; _dirty = true; return *this; }
  Vector4X& x(double x) { _elements[1] = x; _dirty = true; return *this; }
  Vector4X& y(double y) { _elements[2] = y; _dirty = true; return *this; }
  Vector4X& z(double z) { _elements[3] = z; _dirty = true; return *this; }
  double interval() const { return this->invariant(); }
  Vector3X vector3() const { return Vector3X(x(), y(), z()); }
};


class Vector4P : public Vector4 {
public:
  Vector4P(double E, double px, double py, double pz) {
    this->E(E);
    this->px(px);
    this->py(py);
    this->pz(pz);
  }

  ~Vector4P() {}

public:
  double E() const { return _elements[0]; }
  double px() const { return _elements[1]; }
  double py() const { return _elements[2]; }
  double pz() const { return _elements[3]; }
  Vector4P& E(double E)   { _elements[0] = E;  _dirty = true; return *this; }
  Vector4P& px(double px) { _elements[1] = px; _dirty = true; return *this; }
  Vector4P& py(double py) { _elements[2] = py; _dirty = true; return *this; }
  Vector4P& pz(double pz) { _elements[3] = pz; _dirty = true; return *this; }
  double mass() const { return this->invariant(); }
  Vector3P vector3() const { return Vector3P(px(), py(), pz()); }
};


///////////////////////////////////////////////////////

// Define calculations as external functions for flexibility.
// The corresponding member functions are implemented in terms of these.


/// Calculate the modulus-squared of a Vector3,
/// \f$ \sum_{i=1}^3 x_i^2 \f$.
inline double mod2(const Vector3& v3) {
  double mod2 = 0;
  for (size_t i = 0; i < 3; ++i) {
    const double element = v3[i];
    mod2 += element*element;
  }
  return mod2;
}


/// Calculate the modulus of a Vector3,
/// \f$ \sqrt{\sum_{i=1}^3 x_i^2} \f$.
inline double mod(const Vector3& v3) {
  return sqrt(mod2(v3));
}


/// Calculate the Lorentz self-invariant of a Vector4,
/// \f$ v_\mu v^\mu = g_{\mu\nu} x^\mu x^\nu \f$.
inline double invariant(const Vector4& v4) {
  const double metricdiag[4] = {1, -1, -1, -1};
  double invariant = 0;
  for (size_t i = 0; i < 4; ++i) {
    const double element = v4[i];
    invariant += metricdiag[i] * element*element;
  }
  return invariant;
}

inline const string toString(const Vector3& v3) {
  ostringstream out("(");
  out << "("  << v3[0] 
      << ", " << v3[1]
      << ", " << v3[2] << ")";
  return out.str();
}

inline std::ostream& operator<<(std::ostream& out, const Vector3& v3) {
  out << toString(v3);
  return out;
}

inline const string toString(const Vector4& v4) {
  ostringstream out("(");
  out << "("  << v4[0] 
      << "; " << v4[1]
      << ", " << v4[2]
      << ", " << v4[3] << ")";
  return out.str();
}

inline std::ostream& operator<<(std::ostream& out, const Vector4& v4) {
  out << toString(v4);
  return out;
}

