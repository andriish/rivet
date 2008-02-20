// -*- C++ -*-
#ifndef RIVET_MathUtils_HH
#define RIVET_MathUtils_HH

#include "Rivet/Math/MathHeader.hh"
#include <cerrno>

/// A pre-defined value of \f$ \pi \f$.
const double PI = 4*atan(1);
const double TWOPI = 2*PI;

template <typename Real>
inline bool isZero(const Real a, const Real tolerance = 1E-10) {
  return (fabs(a) < tolerance);
}

template <typename Real>
inline bool fuzzyEquals(const Real a, const Real b, const Real tolerance = 1E-10) {
  const double absavg = fabs(a + b)/2.0;
  const double absdiff = fabs(a - b);
  return (absavg == 0.0 && absdiff == 0.0) || absdiff/absavg < tolerance;
}


// Include vectors and matrices
#include "Rivet/Math/Vectors.hh"
#include "Rivet/Math/Matrices.hh"

template <typename Num>
inline Num sqr(Num a) {
  return a*a;
}

// template <typename Num>
// inline Num min(Num a, Num b) {
//   return (a < b) ? a : b;
// }

// template <typename Num>
// inline Num max(Num a, Num b) {
//   return (a < b) ? b : a;
// }

using std::min;
using std::max;

// inline double mapAngle0ToPi(double angle) {}
// inline double mapAngleMPiToPi(double angle) {}
// inline double mapAngle0To2Pi(double angle) {}

inline double delta_phi(double phi1, double phi2) {
  //return min( double(fabs(phi1-phi2)), double(2.*PI-fabs(phi1-phi2)) );
  double deltaPhi = fabs(phi1 - phi2);
  if (deltaPhi > PI) deltaPhi = fabs(deltaPhi - TWOPI);
  assert(deltaPhi >= 0 && deltaPhi <= PI);
  return deltaPhi;
}

inline double delta_rad(double y1, double phi1, double y2, double phi2) {
  const double dphi = delta_phi(phi1, phi2);
  return sqrt( sqr(y1-y2) + sqr(dphi) );
}

inline double deltaR(const Vector3& a, const Vector3& b) {
  return delta_rad(a.pseudorapidity(), a.azimuthalAngle(),
                   b.pseudorapidity(), b.azimuthalAngle());
}

inline double deltaR(const FourVector& a, const FourVector& b, DeltaRScheme scheme = PSEUDORAPIDITY) {
  switch (scheme){
    case PSEUDORAPIDITY :
      return deltaR(a.vector3(), b.vector3());
    case RAPIDITY:
      const FourMomentum* ma = dynamic_cast<const FourMomentum*>(&a);
      const FourMomentum* mb = dynamic_cast<const FourMomentum*>(&b);
      if (!ma || !mb) throw std::runtime_error("deltaR with scheme RAPIDITY, can be called with FourMomenta only");
      return deltaR(*ma, *mb, scheme);
    default: 
      throw std::runtime_error("The specified deltaR scheme is not yet implemented");
  }
}

inline double deltaR(const FourMomentum& a, const FourMomentum& b, DeltaRScheme scheme = PSEUDORAPIDITY) {
  switch (scheme){
    case PSEUDORAPIDITY:
      return deltaR(a.vector3(), b.vector3());
    case RAPIDITY:
      return delta_rad(a.rapidity(), a.azimuthalAngle(),
                       b.rapidity(), b.azimuthalAngle());
    default:
      throw std::runtime_error("The specified deltaR scheme is not yet implemented");
  }
}

// returns phi in the interval (-PI, PI]
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

#endif
