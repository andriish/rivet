// -*- C++ -*-
#ifndef RIVET_MathUtils_HH
#define RIVET_MathUtils_HH

#include "Rivet/Math/MathHeader.hh"
using std::min;
using std::max;

/// @todo Namespace these functions...
//namespace Rivet {

  /// @todo Is the "static" declaration needed?
static const double MAXDOUBLE = std::numeric_limits<double>::max();
static const double MAXINT = std::numeric_limits<int>::max();

  /// A pre-defined value of \f$ \pi \f$.
  const double PI = 4*atan(1);

  /// A pre-defined value of \f$ 2\pi \f$.
  const double TWOPI = 2*PI;

  /// A pre-defined value of \f$ \pi/2 \f$.
  const double HALFPI = PI/2.0;

  /// Enum for signs of numbers.
  enum Sign { MINUS = -1, ZERO = 0, PLUS = 1 };

  /// Compare a floating point number to zero with a degree 
  /// of fuzziness expressed by the absolute @a tolerance parameter.
  inline bool isZero(double val, double tolerance=1E-10) {
    return (fabs(val) < tolerance);
  }

  /// Compare an integral-type number to zero. Since there is no
  /// risk of floating point error, this function just exists in
  /// case @c isZero is accidentally used on an integer type, to avoid 
  /// implicit type conversion. The @a tolerance parameter is ignored.
  inline bool isZero(long val, double tolerance=1E-10) {
    return val == 0;
  }

  /// Compare a number to zero with a degree of fuzziness expressed by the
  /// absolute @a tolerance parameter.
  inline int sign(double val) {
    if (isZero(val)) return ZERO;
    return (val > 0) ? PLUS : MINUS;
  }

  /// Compare a number to zero with a degree of fuzziness expressed by the
  /// absolute @a tolerance parameter.
  inline int sign(long val) {
    if (val == 0) return ZERO;
    return (val > 0) ? PLUS : MINUS;
  }

  /// Compare two floating point numbers with a degree of fuzziness 
  /// expressed by the fractional @a tolerance parameter.
  inline bool fuzzyEquals(double a, double b, double tolerance=1E-10) {
    const double absavg = fabs(a + b)/2.0;
    const double absdiff = fabs(a - b);
    return (absavg == 0.0 && absdiff == 0.0) || absdiff/absavg < tolerance;
  }

  /// Compare two integral-type numbers with a degree of fuzziness.
  /// Since there is no risk of floating point error with integral types, 
  /// this function just exists in case @c fuzzyEquals is accidentally 
  /// used on an integer type, to avoid implicit type conversion. The @a 
  /// tolerance parameter is ignored, even if it would have an
  /// absolute magnitude greater than 1.
  inline bool fuzzyEquals(long a, long b, double tolerance=1E-10) {
    return a == b;
  }

  /// Represents whether an interval is open (non-inclusive) or closed
  /// (inclusive). For example, the interval \f$ [0, \pi) \f$ is closed (an inclusive
  /// boundary) at 0, and open (a non-inclusive boundary) at \f$ \pi \f$.
  enum RangeBoundary { OPEN=0, SOFT=0, CLOSED=1, HARD=1 };

  /// Determine if @a value is in the range @a low to @a high, with boundary
  /// types defined by @a lowbound and @a highbound.
  /// @todo Optimise to one-line at compile time?
  inline bool inRange(double value, double low, double high, 
                      RangeBoundary lowbound=OPEN, RangeBoundary highbound=OPEN) {
    if (lowbound == OPEN && highbound == OPEN) {
      return (value > low && value < high);
    } else if (lowbound == OPEN && highbound == CLOSED) {
      return (value > low && value <= high);
    } else if (lowbound == CLOSED && highbound == OPEN) {
      return (value >= low && value < high);
    } else { // if (lowbound == CLOSED && highbound == CLOSED) {
      return (value >= low && value <= high);
    }
  }
//}


// Include vectors and matrices (can't come earlier?)
#include "Rivet/Math/Vectors.hh"
#include "Rivet/Math/Matrices.hh"


//namespace Rivet {

  /// Named number-type squaring operation.
  template <typename Num>
  inline Num sqr(Num a) {
    return a*a;
  }

  /// Reduce any number to the range [-2PI, 2PI] by repeated addition or
  /// subtraction of 2PI as required. Used to normalise angular measures.
  inline double mapAngleM2PITo2Pi(double angle) {
    return fmod(angle, TWOPI);
  }

  /// Calculate the difference between two angles in radians, returning in the
  /// range [-PI, PI].
  inline double mapAngleMPiToPi(double angle) {
    angle = mapAngleM2PITo2Pi(angle);
    if (angle >  PI) angle -= TWOPI;
    if (angle < -PI) angle += TWOPI;
    assert(angle >= -PI && angle <= PI);
    return angle;
  }

  /// Calculate the difference between two angles in radians, returning in the
  /// range [0, 2PI].
  inline double mapAngle0To2Pi(double angle) {
    angle = mapAngleM2PITo2Pi(angle);
    if (angle < 0) angle += TWOPI;
    assert(angle >= 0 && angle <= TWOPI);
    return angle;
  }

  /// Calculate the difference between two angles in radians, returning in the
  /// range [0, PI].
  inline double mapAngle0ToPi(double angle) {
    return fabs(mapAngleMPiToPi(angle));
  }

  /// Calculate the difference between two angles in radians, returning in the
  /// range [0, PI].
  inline double delta_phi(double phi1, double phi2) {
    return mapAngle0ToPi(fabs(phi1 - phi2));
  }

  /// Calculate the distance between two points in 2D rapidity-azimuthal
  /// ("eta-phi") space. The phi values are given in radians.
  inline double deltaR(double y1, double phi1, double y2, double phi2) {
    const double dphi = delta_phi(phi1, phi2);
    return sqrt( sqr(y1-y2) + sqr(dphi) );
  }

  /// Calculate the distance between two points in 2D rapidity-azimuthal
  /// ("eta-phi") space. The phi values are given in radians.
  inline double delta_rad(double y1, double phi1, double y2, double phi2) {
    return deltaR(y1, phi1, y2, phi2);
  }

  /// Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two
  /// spatial vectors.
  inline double deltaR(const Vector3& a, const Vector3& b) {
    return delta_rad(a.pseudorapidity(), a.azimuthalAngle(),
                     b.pseudorapidity(), b.azimuthalAngle());
  }

  /// Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two
  /// spatial vectors.
  inline double deltaR(const Vector3& v, double eta2, double phi2) {
    return delta_rad(v.pseudorapidity(), v.azimuthalAngle(),
                     eta2, phi2);
  }

  /// Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two
  /// spatial vectors.
  inline double deltaR(double eta1, double phi1, const Vector3& v) {
    return delta_rad(eta1, phi1,
                     v.pseudorapidity(), v.azimuthalAngle());
  }

  /// Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two
  /// four-vectors. There is a scheme ambiguity for momentum-type four vectors
  /// as to whether the pseudorapidity (a purely geometric concept) or the
  /// rapidity (a relativistic energy-momentum quantity) is to be used: this can
  /// be chosen via the optional scheme parameter, which is discouraged in this
  /// case since @c RAPIDITY is only a valid option for vectors whose type is
  /// really the FourMomentum derived class.
  inline double deltaR(const FourVector& a, const FourVector& b, 
                       DeltaRScheme scheme = PSEUDORAPIDITY) {
    switch (scheme) {
    case PSEUDORAPIDITY :
      return deltaR(a.vector3(), b.vector3());
    case RAPIDITY:
      {
        const FourMomentum* ma = dynamic_cast<const FourMomentum*>(&a);
        const FourMomentum* mb = dynamic_cast<const FourMomentum*>(&b);
        if (!ma || !mb) {
          string err = "deltaR with scheme RAPIDITY, can be called with FourMomenta only";
          throw std::runtime_error(err);
        }
        return deltaR(*ma, *mb, scheme);
      }
    default: 
      throw std::runtime_error("The specified deltaR scheme is not yet implemented");
    }
  }


  inline double deltaR(const FourVector& v, 
                       double eta2, double phi2,
                       DeltaRScheme scheme = PSEUDORAPIDITY) {
    switch (scheme) {
    case PSEUDORAPIDITY :
      return deltaR(v.vector3(), eta2, phi2);
    case RAPIDITY:
      {
        const FourMomentum* mv = dynamic_cast<const FourMomentum*>(&v);
        if (!mv) {
          string err = "deltaR with scheme RAPIDITY, can be called with FourMomenta only";
          throw std::runtime_error(err);
        }
        return deltaR(*mv, eta2, phi2, scheme);
      }
    default: 
      throw std::runtime_error("The specified deltaR scheme is not yet implemented");
    }
  }


  inline double deltaR(double eta1, double phi1,
                       const FourVector& v, 
                       DeltaRScheme scheme = PSEUDORAPIDITY) {
    switch (scheme) {
    case PSEUDORAPIDITY :
      return deltaR(eta1, phi1, v.vector3());
    case RAPIDITY:
      {
        const FourMomentum* mv = dynamic_cast<const FourMomentum*>(&v);
        if (!mv) {
          string err = "deltaR with scheme RAPIDITY, can be called with FourMomenta only";
          throw std::runtime_error(err);
        }
        return deltaR(eta1, phi1, *mv, scheme);
      }
    default: 
      throw std::runtime_error("The specified deltaR scheme is not yet implemented");
    }
  }


  /// Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two
  /// four-vectors. There is a scheme ambiguity for momentum-type four vectors
  /// as to whether the pseudorapidity (a purely geometric concept) or the
  /// rapidity (a relativistic energy-momentum quantity) is to be used: this can
  /// be chosen via the optional scheme parameter.
  inline double deltaR(const FourMomentum& a, const FourMomentum& b, 
                       DeltaRScheme scheme = PSEUDORAPIDITY) {
    switch (scheme) {
    case PSEUDORAPIDITY:
      return deltaR(a.vector3(), b.vector3());
    case RAPIDITY:
      return delta_rad(a.rapidity(), a.azimuthalAngle(),
                       b.rapidity(), b.azimuthalAngle());
    default:
      throw std::runtime_error("The specified deltaR scheme is not yet implemented");
    }
  }

  inline double deltaR(const FourMomentum& v,
                       double eta2, double phi2,
                       DeltaRScheme scheme = PSEUDORAPIDITY) {
    switch (scheme) {
    case PSEUDORAPIDITY:
      return deltaR(v.vector3(), eta2, phi2);
    case RAPIDITY:
      return delta_rad(v.rapidity(), v.azimuthalAngle(),
                       eta2, phi2);
    default:
      throw std::runtime_error("The specified deltaR scheme is not yet implemented");
    }
  }


  inline double deltaR(double eta1, double phi1,
                       const FourMomentum& v,
                       DeltaRScheme scheme = PSEUDORAPIDITY) {
    switch (scheme) {
    case PSEUDORAPIDITY:
      return deltaR(eta1, phi1, v.vector3());
    case RAPIDITY:
      return delta_rad(eta1, phi1,
                       v.rapidity(), v.azimuthalAngle());
    default:
      throw std::runtime_error("The specified deltaR scheme is not yet implemented");
    }
  }  


  /// Returns phi in the interval (-PI, PI]
  /// @deprecated Prefer the mapAngleMPiToPi function, whose name makes the operation explicit.
  inline double phi(double px, double py) {
    //return atan2(py, px);
    double value = atan2( py, px );
    if (value > 2*PI || value < -2*PI){
      value = fmod(value, 2*PI);
    }
    if (value <= -PI) value+=2*PI;
    if (value >   PI) value-=2*PI;
    assert(value > -PI && value <= PI);
    return value;
  }


  /// Calculate a rapidity value from the supplied energy @a E and longitudinal momentum @pz.
  inline double rapidity(double E, double pz) {
    if (isZero(E - pz)) {
      throw std::runtime_error("Divergent positive rapidity");
      return MAXDOUBLE;
    }
    if (isZero(E + pz)) {
      throw std::runtime_error("Divergent negative rapidity");
      return -MAXDOUBLE;
    }
    return 0.5*log((E+pz)/(E-pz));
  }


  #include <cerrno>
  /// Calculate a rapidity value from the supplied energy @a E and longitudinal momentum @pz.
  /// @deprecated Taken from D0 and uses error codes unknown to Rivet. 
  /// @todo Doesn't distinguish between positive and negative divergent values.
  inline double y(double E, double pz) {
    errno=0;
    double y;
    if (fabs(E-pz) == 0.) {
      errno=721;
      y = 99999.;
    } else {
      y = 0.5*log((E+pz)/(E-pz));
    }
    return y;
  }

//}

#endif
