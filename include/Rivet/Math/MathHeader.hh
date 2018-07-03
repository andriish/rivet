#ifndef RIVET_Math_MathHeader
#define RIVET_Math_MathHeader

#include "Rivet/Tools/Exceptions.hh"
#include "Rivet/Tools/Utils.hh"

namespace Rivet {


  /// Pre-defined numeric type limits
  /// A pre-defined value of \f$ \pi \f$.
  static const double PI = M_PI;

  /// A pre-defined value of \f$ 2\pi \f$.
  static const double TWOPI = 2*M_PI;

  /// A pre-defined value of \f$ \pi/2 \f$.
  static const double HALFPI = M_PI_2;

  /// Enum for signs of numbers.
  enum Sign { MINUS = -1, ZERO = 0, PLUS = 1 };

  /// Enum for rapidity variable to be used in calculating \f$ R \f$, applying rapidity cuts, etc.
  enum RapScheme { PSEUDORAPIDITY = 0, ETARAP = 0, RAPIDITY = 1, YRAP = 1 };

  /// Enum for range of \f$ \phi \f$ to be mapped into
  enum PhiMapping { MINUSPI_PLUSPI, ZERO_2PI, ZERO_PI };

}

#endif
