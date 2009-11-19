#ifndef RIVET_Math_MathHeader
#define RIVET_Math_MathHeader

#include <stdexcept>
#include <string>
#include <ostream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>

namespace Rivet {

  using std::string;
  using std::ostream;
  using std::ostringstream;
  using std::cout;
  using std::endl;
  using std::pair;
  using std::vector;
  using std::transform;
  using std::min;
  using std::max;

  const double MAXDOUBLE = std::numeric_limits<double>::max();
  const double MAXINT = std::numeric_limits<int>::max();

  /// A pre-defined value of \f$ \pi \f$.
  const double PI = M_PI;

  /// A pre-defined value of \f$ 2\pi \f$.
  const double TWOPI = 2*M_PI;

  /// A pre-defined value of \f$ \pi/2 \f$.
  const double HALFPI = M_PI_2;

  /// Enum for signs of numbers.
  enum Sign { MINUS = -1, ZERO = 0, PLUS = 1 };

  /// Enum for longitudinal variable to be used in calculating \f$ R \f$
  enum DeltaRScheme { PSEUDORAPIDITY, RAPIDITY };

  /// Enum for range of \f$ \phi \f$ to be mapped into
  enum PhiMapping { MINUSPI_PLUSPI, ZERO_2PI };

}

#endif
