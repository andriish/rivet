#ifndef RIVET_MATH_VECTOR3
#define RIVET_MATH_VECTOR3

#include "Rivet/Math/MathConstants.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Math/VectorN.hh"

namespace Rivet {


  class Vector3;
  typedef Vector3 ThreeVector;
  typedef Vector3 V3;

  class Matrix3;

  Vector3 multiply(const double, const Vector3&);
  Vector3 multiply(const Vector3&, const double);
  Vector3 add(const Vector3&, const Vector3&);
  Vector3 operator*(const double, const Vector3&);
  Vector3 operator*(const Vector3&, const double);
  Vector3 operator/(const Vector3&, const double);
  Vector3 operator+(const Vector3&, const Vector3&);
  Vector3 operator-(const Vector3&, const Vector3&);


  /// @brief Three-dimensional specialisation of Vector.
  class Vector3 : public Vector<3> {

    friend class Matrix3;
    friend Vector3 multiply(const double, const Vector3&);
    friend Vector3 multiply(const Vector3&, const double);
    friend Vector3 add(const Vector3&, const Vector3&);
    friend Vector3 subtract(const Vector3&, const Vector3&);

  public:
    Vector3() : Vector<3>() { }

    template<typename V3TYPE>
    Vector3(const V3TYPE& other) {
      this->setX(other.x());
      this->setY(other.y());
      this->setZ(other.z());
    }

    Vector3(const Vector<3>& other) {
      this->setX(other.get(0));
      this->setY(other.get(1));
      this->setZ(other.get(2));
    }

    Vector3(double x, double y, double z) {
      this->setX(x);
      this->setY(y);
      this->setZ(z);
    }

    ~Vector3() { }


  public:

    static Vector3 mkX() { return Vector3(1,0,0); }
    static Vector3 mkY() { return Vector3(0,1,0); }
    static Vector3 mkZ() { return Vector3(0,0,1); }


  public:

    double x() const { return get(0); }
    double y() const { return get(1); }
    double z() const { return get(2); }
    Vector3& setX(double x) { set(0, x); return *this; }
    Vector3& setY(double y) { set(1, y); return *this; }
    Vector3& setZ(double z) { set(2, z); return *this; }


    /// Dot-product with another vector
    double dot(const Vector3& v) const {
      return _vec.dot(v._vec);
    }

    /// Cross-product with another vector
    Vector3 cross(const Vector3& v) const {
      Vector3 result;
      result._vec = _vec.cross(v._vec);
      return result;
    }

    /// Angle in radians to another vector
    double angle(const Vector3& v) const {
      const double localDotOther = unit().dot(v.unit());
      if (localDotOther > 1.0) return 0.0;
      if (localDotOther < -1.0) return M_PI;
      return acos(localDotOther);
    }


    /// Unit-normalized version of this vector.
    Vector3 unitVec() const {
      double md = mod();
      if ( md <= 0.0 ) return Vector3();
      else return *this * 1.0/md;
    }

    /// Synonym for unitVec
    Vector3 unit() const {
      return unitVec();
    }


    /// Polar projection of this vector into the x-y plane
    Vector3 polarVec() const {
      Vector3 rtn = *this;
      rtn.setZ(0.);
      return rtn;
    }
    /// Synonym for polarVec
    Vector3 perpVec() const {
      return polarVec();
    }
    /// Synonym for polarVec
    Vector3 rhoVec() const {
      return polarVec();
    }

    /// Square of the polar radius (
    double polarRadius2() const {
      return x()*x() + y()*y();
    }
    /// Synonym for polarRadius2
    double perp2() const {
      return polarRadius2();
    }
    /// Synonym for polarRadius2
    double rho2() const {
      return polarRadius2();
    }

    /// Polar radius
    double polarRadius() const {
      return sqrt(polarRadius2());
    }
    /// Synonym for polarRadius
    double perp() const {
      return polarRadius();
    }
    /// Synonym for polarRadius
    double rho() const {
      return polarRadius();
    }

    /// @brief Angle subtended by the vector's projection in x-y and the x-axis.
    ///
    /// @note Returns zero in the case of a vector with null x and y components.
    /// @todo Would it be better to return NaN in the null-perp case? Or throw?!
    double azimuthalAngle(const PhiMapping mapping = ZERO_2PI) const {
      // If this has a null perp-vector, return zero rather than let atan2 set an error state
      // This isn't necessary if the implementation supports IEEE floating-point arithmetic (IEC 60559)... are we sure?
      if (x() == 0 && y() == 0) return 0.0; //< Or return nan / throw an exception?
      // Calculate the arctan and return in the requested range
      const double value = atan2( y(), x() );
      return mapAngle(value, mapping);
    }
    /// Synonym for azimuthalAngle.
    double phi(const PhiMapping mapping = ZERO_2PI) const {
      return azimuthalAngle(mapping);
    }

    /// Tangent of the polar angle
    double tanTheta() const {
      return polarRadius()/z();
    }

    /// Angle subtended by the vector and the z-axis.
    double polarAngle() const {
      // Get number beween [0,PI]
      const double polarangle = atan2(polarRadius(), z());
      return mapAngle0ToPi(polarangle);
    }

    /// Synonym for polarAngle
    double theta() const {
      return polarAngle();
    }

    /// @brief Purely geometric approximation to rapidity
    ///
    /// eta = -ln[ tan(theta/2) ]
    ///
    /// Also invariant under z-boosts, equal to y for massless particles.
    ///
    /// Implemented using the tan half-angle formula
    /// tan(theta/2) = sin(theta) / [1 + cos(theta)] = pT / (p + pz)
    double pseudorapidity() const {
      if (mod() == 0.0) return 0.0;
      const double eta = std::log((mod() + fabs(z())) / perp());
      return std::copysign(eta, z());
    }

    /// Synonym for pseudorapidity
    double eta() const {
      return pseudorapidity();
    }

    /// Convenience shortcut for fabs(eta())
    double abseta() const {
      return fabs(eta());
    }


  public:

    /// In-place scalar multiplication operator
    Vector3& operator *= (const double a) {
      _vec = multiply(a, *this)._vec;
      return *this;
    }

    /// In-place scalar division operator
    Vector3& operator /= (const double a) {
      _vec = multiply(1.0/a, *this)._vec;
      return *this;
    }

    /// In-place addition operator
    Vector3& operator += (const Vector3& v) {
      _vec = add(*this, v)._vec;
      return *this;
    }

    /// In-place subtraction operator
    Vector3& operator -= (const Vector3& v) {
      _vec = subtract(*this, v)._vec;
      return *this;
    }

    /// In-place negation operator
    Vector3 operator - () const {
      Vector3 rtn;
      rtn._vec = -_vec;
      return rtn;
    }

  };



  /// Unbound dot-product function
  inline double dot(const Vector3& a, const Vector3& b) {
    return a.dot(b);
  }

  /// Unbound cross-product function
  inline Vector3 cross(const Vector3& a, const Vector3& b) {
    return a.cross(b);
  }

  /// Unbound scalar-product function
  inline Vector3 multiply(const double a, const Vector3& v) {
    Vector3 result;
    result._vec = a * v._vec;
    return result;
  }

  /// Unbound scalar-product function
  inline Vector3 multiply(const Vector3& v, const double a) {
    return multiply(a, v);
  }

  /// Unbound scalar multiplication operator
  inline Vector3 operator * (const double a, const Vector3& v) {
    return multiply(a, v);
  }

  /// Unbound scalar multiplication operator
  inline Vector3 operator * (const Vector3& v, const double a) {
    return multiply(a, v);
  }

  /// Unbound scalar division operator
  inline Vector3 operator / (const Vector3& v, const double a) {
    return multiply(1.0/a, v);
  }

  /// Unbound vector addition function
  inline Vector3 add(const Vector3& a, const Vector3& b) {
    Vector3 result;
    result._vec = a._vec + b._vec;
    return result;
  }

  /// Unbound vector subtraction function
  inline Vector3 subtract(const Vector3& a, const Vector3& b) {
    Vector3 result;
    result._vec = a._vec - b._vec;
    return result;
  }

  /// Unbound vector addition operator
  inline Vector3 operator + (const Vector3& a, const Vector3& b) {
    return add(a, b);
  }

  /// Unbound vector subtraction operator
  inline Vector3 operator - (const Vector3& a, const Vector3& b) {
    return subtract(a, b);
  }

  // More physicsy coordinates etc.

  /// Angle (in radians) between two 3-vectors.
  inline double angle(const Vector3& a, const Vector3& b) {
    return a.angle(b);
  }


  /////////////////////////////////////////////////////


  /// @defgroup momutils_vec3_deta \f$ |\Delta eta| \f$ calculations from 3-vectors
  /// @{

  /// Calculate the difference in pseudorapidity between two spatial vectors.
  inline double deltaEta(const Vector3& a, const Vector3& b, bool sign=false) {
    return deltaEta(a.pseudorapidity(), b.pseudorapidity(), sign);
  }

  /// Calculate the difference in pseudorapidity between two spatial vectors.
  inline double deltaEta(const Vector3& v, double eta2, bool sign=false) {
    return deltaEta(v.pseudorapidity(), eta2, sign);
  }

  /// Calculate the difference in pseudorapidity between two spatial vectors.
  inline double deltaEta(double eta1, const Vector3& v, bool sign=false) {
    return deltaEta(eta1, v.pseudorapidity(), sign);
  }

  /// @}


  /// @defgroup momutils_vec3_dphi \f$ \Delta phi \f$ calculations from 3-vectors
  /// @{

  /// Calculate the difference in azimuthal angle between two spatial vectors.
  inline double deltaPhi(const Vector3& a, const Vector3& b, bool sign=false) {
    return deltaPhi(a.azimuthalAngle(), b.azimuthalAngle(), sign);
  }

  /// Calculate the difference in azimuthal angle between two spatial vectors.
  inline double deltaPhi(const Vector3& v, double phi2, bool sign=false) {
    return deltaPhi(v.azimuthalAngle(), phi2, sign);
  }

  /// Calculate the difference in azimuthal angle between two spatial vectors.
  inline double deltaPhi(double phi1, const Vector3& v, bool sign=false) {
    return deltaPhi(phi1, v.azimuthalAngle(), sign);
  }

  /// @}


  /// @defgroup momutils_vec3_dr \f$ \Delta R \f$ calculations from 3-vectors
  /// @{

  /// Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two spatial vectors.
  inline double deltaR2(const Vector3& a, const Vector3& b) {
    return deltaR2(a.pseudorapidity(), a.azimuthalAngle(),
                   b.pseudorapidity(), b.azimuthalAngle());
  }

  /// Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two spatial vectors.
  inline double deltaR(const Vector3& a, const Vector3& b) {
    return sqrt(deltaR2(a,b));
  }

  /// Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two spatial vectors.
  inline double deltaR2(const Vector3& v, double eta2, double phi2) {
    return deltaR2(v.pseudorapidity(), v.azimuthalAngle(), eta2, phi2);
  }

  /// Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two spatial vectors.
  inline double deltaR(const Vector3& v, double eta2, double phi2) {
    return sqrt(deltaR2(v, eta2, phi2));
  }

  /// Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two spatial vectors.
  inline double deltaR2(double eta1, double phi1, const Vector3& v) {
    return deltaR2(eta1, phi1, v.pseudorapidity(), v.azimuthalAngle());
  }

  /// Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two spatial vectors.
  inline double deltaR(double eta1, double phi1, const Vector3& v) {
    return sqrt(deltaR2(eta1, phi1, v));
  }

  /// @}


  /// @defgroup momutils_vec3_mt MT calculation
  /// @{

  /// Calculate transverse mass of a visible and an invisible 3-vector
  inline double mT(const Vector3& vis, const Vector3& invis) {
    // return sqrt(2*vis.perp()*invis.perp() * (1 - cos(deltaPhi(vis, invis))) );
    return mT(vis.perp(), invis.perp(), deltaPhi(vis, invis));
  }

  /// @}


}

#endif
