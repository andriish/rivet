#ifndef MATH_MATHHEADER
#define MATH_MATHHEADER

#include <stdexcept>
#include <string>
#include <ostream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <map>
#include <vector>

#include "Rivet/Math/eigen/vector.h"
#include "Rivet/Math/eigen/matrix.h"

using std::string;
using std::ostream;
using std::ostringstream;
using std::cout;
using std::endl;
using std::pair;
using std::vector;

enum DeltaRScheme {PSEUDORAPIDITY, RAPIDITY};
enum PhiMapping {MINUSPIPLUSPI, ZERO2PI};

#endif
