#ifndef RIVET_Rivet_H
#define RIVET_Rivet_H

#include <typeinfo>
#include <set>
#include <list>
#include <map>
#include <utility>
#include <string>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <cassert>

/// This is the main namespace in which all Rivet classes are defined.
namespace Rivet {

  // Convenient imports of common STL classes and functions.
  using std::set;
  using std::map;
  using std::multimap;
  using std::type_info;
  using std::string;
  using std::stringstream;
  using std::less;
  using std::list;
  using std::vector;
  using std::pair;
  using std::make_pair;
  using std::runtime_error;
  using std::min;
  using std::max;
  using std::numeric_limits;
  using std::ostream;
  using std::istream;
  using std::cout;
  using std::cin;
  using std::cerr;
  using std::setw;
  using std::endl;

  /// A sensible default maximum value of rapidity for Rivet analyses to use.
  const double MaxRapidity = 100000.0;

  /// A pre-defined value of pi.
  const double PI = 4*atan(1);

  /// Convenient function for streaming out vectors of any streamable object.
  template<typename T>
  inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    os << "[ ";
    for (typename std::vector<T>::const_iterator i = vec.begin(); i != vec.end(); ++i) {
      os << *i << " ";
    }
    os << "]";
    return os;
  }

}

#endif
