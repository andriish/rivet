#ifndef MATH_VECTOR4
#define MATH_VECTOR4

#include "Rivet/Math/MathHeader.hh"
#include "Rivet/Math/VectorN.hh"
#include "Rivet/Math/Vector3.hh"

class FourVector;
class FourMomentum;
class LorentzTransform;
typedef FourVector Vector4;
const double invariant(const FourVector& lv);
FourVector transform(const LorentzTransform& lt, const FourVector& v4);


class FourVector : public Vector<4> {
  friend FourVector multiply(const double a, const FourVector& v);
  friend FourVector multiply(const FourVector& v, const double a);
  friend FourVector add(const FourVector& a, const FourVector& b);
  friend FourVector transform(const LorentzTransform& lt, const FourVector& v4);

public:
  FourVector() { }

  template<typename V4>
  FourVector(const V4& other) {
    this->t(other.t());
    this->x(other.x());
    this->y(other.y());
    this->z(other.z());
  }

  FourVector(const Vector<4>& other) 
  : Vector<4>(other) { }

  FourVector(const double t, const double x, const double y, const double z) {
    this->t(t);
    this->x(x);
    this->y(y);
    this->z(z);
  }

  virtual ~FourVector() { }

public:
  const double t() const { return get(0); }
  const double x() const { return get(1); }
  const double y() const { return get(2); }
  const double z() const { return get(3); }
  FourVector& t(const double t) { set(0, t); return *this; }
  FourVector& x(const double x) { set(1, x); return *this; }
  FourVector& y(const double y) { set(2, y); return *this; }
  FourVector& z(const double z) { set(3, z); return *this; }
  
  double invariant() const {
    return ::invariant(*this);
  }

  double angle(const FourVector& v) const {
    return ::angle(*this, v);
  }

  double angle(const Vector3& v3) const {
    return ::angle(*this, v3);
  }

  double polarRadius2() const {
    return ::polarRadius2(*this);
  }

  double polarRadius() const {
    return ::polarRadius(*this);    
  }
  
  double azimuthalAngle() const {
    return ::azimuthalAngle(*this);
  }
  
  double polarAngle() const {
    return ::polarAngle(*this);
  }
  
  double pseudorapidity() const {
    return ::pseudorapidity(*this);
  }

  /// Get the spatial part of the 4-vector as a 3-vector.
  Vector3 vector3() const {
    return Vector3(get(1), get(2), get(3));
  }

public:
  FourVector& operator*=(double a) {
    _vec = multiply(a, *this)._vec;
    return *this;
  }

  FourVector& operator/=(double a) {
    _vec = multiply(1.0/a, *this)._vec;
    return *this;
  }

  FourVector& operator+=(FourVector v) {
    _vec = add(*this, v)._vec;
    return *this;
  }

  FourVector& operator-=(FourVector v) {
    _vec = add(*this, -v)._vec;
    return *this;
  }

  FourVector operator-() const {
    FourVector result;
    result._vec = -_vec;
    return result;
  }

};


/// Contract two 4-vectors, with metric signature (+ - - -).
inline double contract(const FourVector& a, const FourVector& b) {
  double result = a.t()*b.t() - a.x()*b.x() - a.y()*b.y() - a.z()*b.z();
  return result;
}

/// Contract two 4-vectors, with metric signature (+ - - -).
inline double dot(const FourVector& a, const FourVector& b) {
  return contract(a, b);
}

/// Contract two 4-vectors, with metric signature (+ - - -).
inline double operator*(const FourVector& a, const FourVector& b) {
  return contract(a, b);
}

inline FourVector multiply(const double a, const FourVector& v) {
  FourVector result;
  result._vec = a * v._vec;
  return result;
}

inline FourVector multiply(const FourVector& v, const double a) {
  return multiply(a, v);
}

inline FourVector operator*(const double a, const FourVector& v) {
  return multiply(a, v);
}

inline FourVector operator*(const FourVector& v, const double a) {
  return multiply(a, v);
}

inline FourVector operator/(const FourVector& v, const double a) {
  return multiply(1.0/a, v);
}

inline FourVector add(const FourVector& a, const FourVector& b) {
  FourVector result;
  result._vec = a._vec + b._vec;
  return result;
}

inline FourVector operator+(const FourVector& a, const FourVector& b) {
  return add(a, b);
}

inline FourVector operator-(const FourVector& a, const FourVector& b) {
  return add(a, -b);
}

/// Calculate the Lorentz self-invariant of a 4-vector.
/// \f$ v_\mu v^\mu = g_{\mu\nu} x^\mu x^\nu \f$.
inline const double invariant(const FourVector& lv) {
  const double inv = lv.t()*lv.t() - lv.x()*lv.x() - lv.y()*lv.y() - lv.z()*lv.z();
  return inv;
}

/// Angle (in radians) between spatial parts of two Lorentz vectors.
inline double angle(const FourVector& a, const FourVector& b) {
  return angle( a.vector3(), b.vector3() );
}

/// Angle (in radians) between spatial parts of two Lorentz vectors.
inline double angle(const Vector3& a, const FourVector& b) {
  return angle( a, b.vector3() );
}

/// Angle (in radians) between spatial parts of two Lorentz vectors.
inline double angle(const FourVector& a, const Vector3& b) {
  return angle( a.vector3(), b );
}

/// Calculate transverse length sq. \f$ \rho^2 \f$ of a Lorentz vector.
inline double polarRadius2(const FourVector& v) {
  return polarRadius2(v.vector3());
}
    
/// Calculate transverse length \f$ \rho \f$ of a Lorentz vector.
inline double polarRadius(const FourVector& v) {
  return polarRadius(v.vector3());
}

/// Calculate azimuthal angle of a Lorentz vector.
inline double azimuthalAngle(const FourVector& v) {
  return azimuthalAngle(v.vector3());
}

/// Calculate polar angle of a Lorentz vector.
inline double polarAngle(const FourVector& v) {
  return polarAngle(v.vector3());
}

/// Calculate pseudorapidity of a Lorentz vector.
inline double pseudorapidity(const FourVector& v) {
  return pseudorapidity(v.vector3());
}



////////////////////////////////////////////////


double rapidity(const FourMomentum& v);
double mass2(const FourMomentum& v);
double mass(const FourMomentum& v);
double pT2(const FourMomentum& lv);
double pT(const FourMomentum& lv);
double Et(const FourMomentum& lv);


/// Specialized version of the FourVector with momentum/energy functionality.
class FourMomentum : public FourVector {
public:
  FourMomentum() { }

  template<typename V4>
  FourMomentum(const V4& other) {
    this->E(other.t());
    this->px(other.x());
    this->py(other.y());
    this->pz(other.z());
  }

  FourMomentum(const Vector<4>& other) 
    : FourVector(other) { }
  
  FourMomentum(const double E, const double px, const double py, const double pz) {
    this->E(E);
    this->px(px);
    this->py(py);
    this->pz(pz);
  }

  ~FourMomentum() {}

public:
  /// Get energy \f$ E \f$ (time component of momentum).
  double E() const { return t(); }

  /// Get 3-momentum part, \f$ p \f$.
  Vector3 p() const { return vector3(); }

  /// Get x-component of momentum \f$ p_x \f$.
  double px() const { return x(); }

  /// Get y-component of momentum \f$ p_y \f$.
  double py() const { return y(); }

  /// Get z-component of momentum \f$ p_z \f$.
  double pz() const { return z(); }

  /// Set energy \f$ E \f$ (time component of momentum).
  FourMomentum& E(double E)   { t(E); return *this; }

  /// Set x-component of momentum \f$ p_x \f$.
  FourMomentum& px(double px) { x(px); return *this; }

  /// Set y-component of momentum \f$ p_y \f$.
  FourMomentum& py(double py) { y(py); return *this; }

  /// Set z-component of momentum \f$ p_z \f$.
  FourMomentum& pz(double pz) { z(pz); return *this; }

  /// Get squared mass \f$ m^2 = E^2 - p^2 \f$ (the Lorentz self-invariant).
  double mass2() const { return ::mass2(*this); }

  /// Get mass \f$ m = \sqrt{E^2 - p^2} \f$ (the Lorentz self-invariant).
  double mass() const { return ::mass(*this); }

  /// Calculate rapidity.
  double rapidity() const { return ::rapidity(*this); }

  /// Calculate squared transverse momentum \f$ p_T^2 \f$.
  double pT2() const { return ::pT2(*this); }

  /// Calculate transverse momentum \f$ p_T \f$.
  double pT() const { return ::pT(*this); }

  /// Calculate transverse energy \f$ E_T = E \cos{\theta} \f$.
  double Et() const { return ::Et(*this); }

  /// Calculate boost vector.
  Vector3 boostVector() const { return Vector3(px()/E(), py()/E(), pz()/E()); }
};


/// Get squared mass \f$ m^2 = E^2 - p^2 \f$ (the Lorentz self-invariant) of a momentum 4-vector.
inline double mass2(const FourMomentum& v) { 
  return invariant(v); 
}

/// Get mass \f$ m = \sqrt{E^2 - p^2} \f$ (the Lorentz self-invariant) of a momentum 4-vector.
inline double mass(const FourMomentum& v) { 
  assert(mass2(v) >= 0);
  return sqrt(mass2(v)); 
}

/// Calculate rapidity of a momentum 4-vector.
inline double rapidity(const FourMomentum& v) {
  return 0.5 * log( (v.E() + v.pz()) / (v.E() - v.pz()) );
}

/// Calculate squared transverse momentum \f$ p_T^2 \f$ of a momentum 4-vector.
inline double pT2(const FourMomentum& lv) {
  return polarRadius2(lv.vector3());
}
    
/// Calculate transverse momentum \f$ p_T \f$ of a momentum 4-vector.
inline double pT(const FourMomentum& lv) {
  return sqrt(pT2(lv));
}

/// Calculate transverse energy \f$ E_T = E \cos{\theta} \f$ of a momentum 4-vector.
inline double Et(const FourMomentum& lv) {
  return sin(lv.polarAngle()) * lv.E();
}

/// Calculate velocity boost vector of a momentum 4-vector.
inline Vector3 boostVector(const FourMomentum& lv) {
  const Vector3 p3 = lv.vector3();
  const double m2 = lv.mass2();
  if (isZero(m2)) return p3.unit();
  else {
    // Could also do this via beta = tanh(rapidity), but that's
    // probably more messy from a numerical hygiene point of view.
    const double p2 = p3.mod2();
    const double beta = sqrt( p2 / (m2 + p2) );
    return beta * p3.unit();
  }
}


//////////////////////////////////////////////////////


/// Render a 4-vector as a string.
inline const string toString(const FourVector& lv) {
  ostringstream out;
  out << "("  << (fabs(lv.t()) < 1E-30 ? 0.0 : lv.t())
      << "; " << (fabs(lv.x()) < 1E-30 ? 0.0 : lv.x())
      << ", " << (fabs(lv.y()) < 1E-30 ? 0.0 : lv.y())
      << ", " << (fabs(lv.z()) < 1E-30 ? 0.0 : lv.z())
      << ")";
  return out.str();
}

/// Write a 4-vector to an ostream.
inline std::ostream& operator<<(std::ostream& out, const FourVector& lv) {
  out << toString(lv);
  return out;
}

#endif
