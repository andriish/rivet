#ifndef RIVET_Rivet_H
#define RIVET_Rivet_H

#include <typeinfo>
#include <set>
#include <list>
#include <map>
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

  const double MaxRapidity = 100000.0;
  const double PI = 4*atan(1);

}

#endif
