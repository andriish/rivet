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

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"


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

  using HepMC::GenEvent;
  using HepMC::GenParticle;
  using HepMC::GenVertex;

  /// A sensible default maximum value of rapidity for Rivet analyses to use.
  const double MaxRapidity = 100000.0;

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

// Now import some Rivet classes
#include "Rivet/RivetAIDA.fhh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Math/Vectors.hh"
#include "Rivet/Math/Matrices.hh"
#include "Rivet/Math/Units.hh"

#endif
