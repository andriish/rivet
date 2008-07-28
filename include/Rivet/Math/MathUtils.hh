// -*- C++ -*-
#ifndef RIVET_MathUtils_HH
#define RIVET_MathUtils_HH

#include "Rivet/Math/MathHeader.hh"

namespace Rivet {


  /// Compare a floating point number to zero with a degree 
  /// of fuzziness expressed by the absolute @a tolerance parameter.
  inline bool isZero(double val, double tolerance=1E-8) {
    return (fabs(val) < tolerance);
  }

  /// Compare an integral-type number to zero. Since there is no
  /// risk of floating point error, this function just exists in
  /// case @c isZero is accidentally used on an integer type, to avoid 
  /// implicit type conversion. The @a tolerance parameter is ignored.
  inline bool isZero(long val, double tolerance=1E-8) {
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
  inline bool fuzzyEquals(double a, double b, double tolerance=1E-5) {
    const double absavg = fabs(a + b)/2.0;
    const double absdiff = fabs(a - b);
    const bool rtn = (absavg == 0.0 && absdiff == 0.0) || absdiff/absavg < tolerance;
    return rtn;
  }

  /// Compare two integral-type numbers with a degree of fuzziness.
  /// Since there is no risk of floating point error with integral types, 
  /// this function just exists in case @c fuzzyEquals is accidentally 
  /// used on an integer type, to avoid implicit type conversion. The @a 
  /// tolerance parameter is ignored, even if it would have an
  /// absolute magnitude greater than 1.
  inline bool fuzzyEquals(long a, long b, double tolerance=1E-5) {
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

  /// Named number-type squaring operation.
  template <typename Num>
  inline Num sqr(Num a) {
    return a*a;
  }

  /// Reduce any number to the range [-2PI, 2PI] by repeated addition or
  /// subtraction of 2PI as required. Used to normalise angular measures.
  inline double _mapAngleM2PITo2Pi(double angle) {
    double rtn = fmod(angle, TWOPI);
    assert(rtn >= -TWOPI && rtn <= TWOPI);
    return rtn;
  }

  /// Calculate the difference between two angles in radians, returning in the
  /// range (-PI, PI].
  inline double mapAngleMPiToPi(double angle) {
    double rtn = _mapAngleM2PITo2Pi(angle);
    rtn = (rtn >   PI ? rtn-TWOPI :
           rtn <= -PI ? rtn+TWOPI : rtn);
    assert(rtn > -PI && rtn <= PI);
    return rtn;
  }

  /// Calculate the difference between two angles in radians, returning in the
  /// range [0, 2PI).
  inline double mapAngle0To2Pi(double angle) {
    double rtn = _mapAngleM2PITo2Pi(angle);
    if (rtn < 0) rtn += TWOPI;
    if (rtn == TWOPI) rtn = 0;
    assert(rtn >= 0 && rtn < TWOPI);
    return rtn;
  }

  /// Calculate the difference between two angles in radians, returning in the
  /// range [0, PI].
  inline double mapAngle0ToPi(double angle) {
    double rtn = fabs(mapAngleMPiToPi(angle));
    assert(rtn >= 0 && rtn <= PI);
    return rtn;
  }

  /// Calculate the difference between two angles in radians, 
  /// returning in the range [0, PI].
  /// @todo Not suitable for calculating polar angles
  inline double deltaPhi(double phi1, double phi2) {
    double angle = fabs(mapAngleMPiToPi(phi1 - phi2));
    assert(angle >= 0 && angle <= PI);
    return angle;
  }

  /// Calculate the distance between two points in 2D rapidity-azimuthal
  /// ("eta-phi") space. The phi values are given in radians.
  inline double deltaR(double y1, double phi1, double y2, double phi2) {
    const double dphi = deltaPhi(phi1, phi2);
    return sqrt( sqr(y1-y2) + sqr(dphi) );
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


  /// Calculate (pseudo)rapidity relative to another vector
  /// @todo Implement in Vector{3,4}
  //inline double rapidity(v4, relativeTo) { ... }


  //////////////////////////////////////////////////////////////
  /// DEPRECATED
  //////////////////////////////////////////////////////////////

  /// Calculate the difference between two angles in radians, 
  /// returning in the range [0, PI].
  /// @deprecated Prefer deltaPhi
  inline double delta_phi(double phi1, double phi2) {
    return deltaPhi(phi1, phi2);
  }


  /// Calculate the distance between two points in 2D rapidity-azimuthal
  /// ("eta-phi") space. The phi values are given in radians.
  /// @deprecated Prefer deltaR(y1, phi1, y2, phi2)
  inline double delta_rad(double y1, double phi1, double y2, double phi2) {
    return deltaR(y1, phi1, y2, phi2);
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


}

#endif
