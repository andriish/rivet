#ifndef MATH_VECTOR3
#define MATH_VECTOR3

#include "Rivet/Math/MathHeader.hh"
#include "Rivet/Math/VectorN.hh"

class Vector3;
typedef Vector3 ThreeVector;
class Matrix3;

double dot(const Vector3&, const Vector3&);
Vector3 cross(const Vector3&, const Vector3&);
Vector3 multiply(const double, const Vector3&);
Vector3 multiply(const Vector3&, const double);
Vector3 add(const Vector3&, const Vector3&);
Vector3 operator*(const double, const Vector3&);
Vector3 operator*(const Vector3&, const double);
Vector3 operator/(const Vector3&, const double);
Vector3 operator+(const Vector3&, const Vector3&);
Vector3 operator-(const Vector3&, const Vector3&);
double angle(const Vector3&, const Vector3&);
double polarRadius2(const Vector3&);
double polarRadius(const Vector3&);
double azimuthalAngle(const Vector3&, const PhiMapping mapping = MINUSPIPLUSPI);
double polarAngle(const Vector3&);
double pseudorapidity(const Vector3&);


class Vector3 : public Vector<3> {
  friend class Matrix3;
  friend double dot(const Vector3&, const Vector3&);
  friend Vector3 cross(const Vector3&, const Vector3&);
  friend Vector3 multiply(const double, const Vector3&);
  friend Vector3 multiply(const Vector3&, const double);
  friend Vector3 add(const Vector3&, const Vector3&);
  friend Vector3 subtract(const Vector3&, const Vector3&);

public:
  Vector3() : Vector<3>() { }

  template<typename V3>
  Vector3(const V3& other) {
    this->x(other.x());
    this->y(other.y());
    this->z(other.z());
  }

  Vector3(const Vector<3>& other) {
    this->x(other.get(0));
    this->y(other.get(1));
    this->z(other.get(2));
  }

  Vector3(double x, double y, double z) {
    this->x(x);
    this->y(y);
    this->z(z);
  }

  ~Vector3() { }

public:
  static Vector3 X() { return Vector3(1,0,0); }
  static Vector3 Y() { return Vector3(0,1,0); }
  static Vector3 Z() { return Vector3(0,0,1); }

public:
  double x() const { return get(0); }
  double y() const { return get(1); }
  double z() const { return get(2); }
  Vector3& x(const double x) { set(0, x); return *this; }
  Vector3& y(const double y) { set(1, y); return *this; }
  Vector3& z(const double z) { set(2, z); return *this; }

  double mod2() const {
    return ::mod2(*this);
  }

  double mod() const { 
    return ::mod(*this);
  }

  double dot(const Vector3& v) const { 
    return ::dot(*this, v);
  }

  Vector3 cross(const Vector3& v) const { 
    return ::cross(*this, v);
  }

  double angle(const Vector3& v) const {
    return ::angle(*this, v);
  }

  Vector3 unit() const {
    if (isZero()) return *this;
    else return *this * 1.0/this->mod();
  }

  double polarRadius2() const {
    return ::polarRadius2(*this);
  }

  double polarRadius() const {
    return ::polarRadius(*this);    
  }
  
  double azimuthalAngle(const PhiMapping mapping = MINUSPIPLUSPI) const {
    return ::azimuthalAngle(*this, mapping);
  }
  
  double polarAngle() {
    return ::polarAngle(*this);
  }
  
  double pseudorapidity() const {
    return ::pseudorapidity(*this);
  }


public:
  Vector3& operator*=(const double a) {
    _vec = multiply(a, *this)._vec;
    return *this;
  }

  Vector3& operator/=(const double a) {
    _vec = multiply(1.0/a, *this)._vec;
    return *this;
  }

  Vector3& operator+=(const Vector3& v) {
    _vec = add(*this, v)._vec;
    return *this;
  }

  Vector3& operator-=(const Vector3& v) {
    _vec = subtract(*this, v)._vec;
    return *this;
  }

  Vector3 operator-() const {
    Vector3 rtn;
    rtn._vec = -_vec;
    return rtn;
  }

};



inline double dot(const Vector3& a, const Vector3& b) {
  return a._vec.dot(b._vec);
}

inline Vector3 cross(const Vector3& a, const Vector3& b) {
  Vector3 result;
  result._vec = a._vec.cross(b._vec);
  //cout << "&&&" << a << " x " << b << " -> " << result._vec << endl;
  return result;
}

inline Vector3 multiply(const double a, const Vector3& v) {
  Vector3 result;
  result._vec = a * v._vec;
  return result;
}

inline Vector3 multiply(const Vector3& v, const double a) {
  return multiply(a, v);
}

inline Vector3 operator*(const double a, const Vector3& v) {
  return multiply(a, v);
}

inline Vector3 operator*(const Vector3& v, const double a) {
  return multiply(a, v);
}

inline Vector3 operator/(const Vector3& v, const double a) {
  return multiply(1.0/a, v);
}

inline Vector3 add(const Vector3& a, const Vector3& b) {
  Vector3 result;
  result._vec = a._vec + b._vec;
  return result;
}

inline Vector3 subtract(const Vector3& a, const Vector3& b) {
  Vector3 result;
  result._vec = a._vec - b._vec;
  return result;
}

inline Vector3 operator+(const Vector3& a, const Vector3& b) {
  return add(a, b);
}

inline Vector3 operator-(const Vector3& a, const Vector3& b) {
  return subtract(a, b);
}

/// More physicsy coordinates etc.

/// Angle (in radians) between two 3-vectors.
inline double angle(const Vector3& a, const Vector3& b) {
  return acos( dot(a.unit(), b.unit()) );
}

/// Calculate transverse length sq. \f$ \rho^2 \f$ of a 3-vector.
inline double polarRadius2(const Vector3& v) {
  return v.x()*v.x() + v.y()*v.y();
}
    
/// Calculate transverse length \f$ \rho \f$ of a 3-vector.
inline double polarRadius(const Vector3& v) {
  return sqrt(polarRadius2(v));
}

/// Calculate azimuthal angle of a 3-vector.
/// returns a number in (-pi, pi]
/// or in [0, 2pi) according to the mapping scheme selected
inline double azimuthalAngle(const Vector3& v, const PhiMapping mapping) {
  //return atan( v.y() / v.x() );
  double value = atan2( v.y(), v.x() );
  if (value > 2*PI || value < -2*PI){
    value = fmod(value, 2*PI);
  }
  if (value <= -PI) value+=2*PI;
  if (value >   PI) value-=2*PI;
  
  switch (mapping){
    case MINUSPIPLUSPI:
      assert(value > -PI && value <= PI);
      return value;
    case ZERO2PI:
      if (value >= 0) {
        assert(value >= 0 && value < 2*PI);
        return value;
      }
      else {
        value = 2*PI + value;
        assert(value >= 0 && value < 2*PI);
        return value;
      }
    default:
      throw std::runtime_error("The specified Phi mapping scheme is not yet implemented"); 
  }
}

/// Calculate polar angle of a 3-vector.
inline double polarAngle(const Vector3& v) {
  double polarangle = atan( polarRadius(v) / v.z() );
  if (polarangle < 0) polarangle += PI;
  return polarangle;
}

/// Calculate pseudorapidity of a 3-vector.
inline double pseudorapidity(const Vector3& v) {
  return -log(tan( 0.5 * polarAngle(v) ));
}

#endif
