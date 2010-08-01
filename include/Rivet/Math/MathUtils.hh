// -*- C++ -*-
#ifndef RIVET_MathUtils_HH
#define RIVET_MathUtils_HH

#include "Rivet/Math/MathHeader.hh"
#include <cassert>

namespace Rivet {


  /// @name Number comparisons etc.
  //@{

  /// Compare a floating point number to zero with a degree
  /// of fuzziness expressed by the absolute @a tolerance parameter.
  inline bool isZero(double val, double tolerance=1E-8) {
    return (fabs(val) < tolerance);
  }

  /// Compare an integral-type number to zero. Since there is no
  /// risk of floating point error, this function just exists in
  /// case @c isZero is accidentally used on an integer type, to avoid
  /// implicit type conversion. The @a tolerance parameter is ignored.
  inline bool isZero(long val, double UNUSED(tolerance)=1E-8) {
    return val == 0;
  }

  /// Find the sign of a number
  inline int sign(double val) {
    if (isZero(val)) return ZERO;
    const int valsign = (val > 0) ? PLUS : MINUS;
    return valsign;
  }

  /// Find the sign of a number
  inline int sign(int val) {
    if (val == 0) return ZERO;
    return (val > 0) ? PLUS : MINUS;
  }

  /// Find the sign of a number
  inline int sign(long val) {
    if (val == 0) return ZERO;
    return (val > 0) ? PLUS : MINUS;
  }

  /// Compare two floating point numbers with a degree of fuzziness
  /// expressed by the fractional @a tolerance parameter.
  inline bool fuzzyEquals(double a, double b, double tolerance=1E-5) {
    const double absavg = fabs(a + b)/2.0;
    const double absdiff = fabs(a - b);
    const bool rtn = (absavg == 0.0 && absdiff == 0.0) || absdiff < tolerance*absavg;
    return rtn;
  }

  /// Compare two integral-type numbers with a degree of fuzziness.
  /// Since there is no risk of floating point error with integral types,
  /// this function just exists in case @c fuzzyEquals is accidentally
  /// used on an integer type, to avoid implicit type conversion. The @a
  /// tolerance parameter is ignored, even if it would have an
  /// absolute magnitude greater than 1.
  inline bool fuzzyEquals(long a, long b, double UNUSED(tolerance)=1E-5) {
    return a == b;
  }

  /// Represents whether an interval is open (non-inclusive) or closed
  /// (inclusive). For example, the interval \f$ [0, \pi) \f$ is closed (an inclusive
  /// boundary) at 0, and open (a non-inclusive boundary) at \f$ \pi \f$.
  enum RangeBoundary { OPEN=0, SOFT=0, CLOSED=1, HARD=1 };

  /// Determine if @a value is in the range @a low to @a high, with boundary
  /// types defined by @a lowbound and @a highbound.
  /// @todo Optimise to one-line at compile time?
  template<typename NUM>
  inline bool inRange(NUM value, NUM low, NUM high,
                      RangeBoundary lowbound=CLOSED, RangeBoundary highbound=OPEN) {
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


  /// Determine if @a value is in the range @a low to @a high, with boundary
  /// types defined by @a lowbound and @a highbound.
  /// @todo Optimise to one-line at compile time?
  inline bool inRange(int value, int low, int high,
                      RangeBoundary lowbound=CLOSED, RangeBoundary highbound=CLOSED) {
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

  //@}



  /// @name Statistics functions
  //@{

  /// Calculate the mean of a sample
  inline double mean(const vector<int>& sample) {
    double mean = 0.0;
    for (size_t i=0; i<sample.size(); ++i) {
      mean += sample[i];
    }
    return mean/sample.size();
  }


  /// Calculate the covariance (variance) between two samples
  inline double covariance(const vector<int>& sample1, const vector<int>& sample2) {
    double mean1 = mean(sample1);
    double mean2 = mean(sample2);
    int N = sample1.size();
    double cov = 0.0;
    for (int i = 0; i < N; i++) {
      double cov_i = (sample1[i] - mean1)*(sample2[i] - mean2);
      cov += cov_i;
    }
    if (N > 1) return cov/(N-1);
    else return 0.0;
  }


  /// Calculate the correlation strength between two samples
  inline double correlation(const vector<int>& sample1, const vector<int>& sample2) {
    const double cov = covariance(sample1, sample2);
    const double var1 = covariance(sample1, sample1);
    const double var2 = covariance(sample2, sample2);
    const double correlation = cov/sqrt(var1*var2);
    const double corr_strength = correlation*sqrt(var2/var1);
    return corr_strength;
  }

  //@}


  /// @name Angle range mappings
  //@{

  /// Reduce any number to the range [-2PI, 2PI] by repeated addition or
  /// subtraction of 2PI as required. Used to normalise angular measures.
  inline double _mapAngleM2PITo2Pi(double angle) {
    double rtn = fmod(angle, TWOPI);
    if (isZero(rtn)) return 0;
    assert(rtn >= -TWOPI && rtn <= TWOPI);
    return rtn;
  }

  /// Map an angle into the range (-PI, PI].
  inline double mapAngleMPiToPi(double angle) {
    double rtn = _mapAngleM2PITo2Pi(angle);
    if (isZero(rtn)) return 0;
    rtn = (rtn >   PI ? rtn-TWOPI :
           rtn <= -PI ? rtn+TWOPI : rtn);
    assert(rtn > -PI && rtn <= PI);
    return rtn;
  }

  /// Map an angle into the range [0, 2PI).
  inline double mapAngle0To2Pi(double angle) {
    double rtn = _mapAngleM2PITo2Pi(angle);
    if (isZero(rtn)) return 0;
    if (rtn < 0) rtn += TWOPI;
    if (rtn == TWOPI) rtn = 0;
    assert(rtn >= 0 && rtn < TWOPI);
    return rtn;
  }

  /// Map an angle into the range [0, PI].
  inline double mapAngle0ToPi(double angle) {
    double rtn = fabs(mapAngleMPiToPi(angle));
    if (isZero(rtn)) return 0;
    assert(rtn > 0 && rtn <= PI);
    return rtn;
  }

  //@}


  /// @name Phase space measure helpers
  //@{

  /// Calculate the difference between two angles in radians,
  /// returning in the range [0, PI].
  inline double deltaPhi(double phi1, double phi2) {
    return mapAngle0ToPi(phi1 - phi2);
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


}

//@}

#endif
