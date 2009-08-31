#ifndef RIVET_Rivet_HH
#define RIVET_Rivet_HH

#include "Rivet/Config/RivetConfig.hh"
#include "Rivet/Config/BuildOptions.hh"

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
#include <fstream>

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
  static const double MAXRAPIDITY = 100000.0;

  /// A function to get the Rivet version string
  string version();

}


// AIDA headers
#include "Rivet/RivetAIDA.fhh"

// Pull some Boost defns into the Rivet namespace
#include "Rivet/RivetBoost.hh"

// HepMC headers and helper functions
#include "Rivet/RivetHepMC.hh"

// Now import some Rivet classes
#include "Rivet/Exceptions.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Math/Vectors.hh"
#include "Rivet/Math/Matrices.hh"
#include "Rivet/Math/Units.hh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/Particle.hh"
#include "Rivet/ParticleName.hh"


namespace Rivet {


  /// Convenient function for streaming out vectors of any streamable object.
  template<typename T>
  inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    os << "[ ";
    foreach (const T& i, vec) {
      os << i << " ";
    }
    os << "]";
    return os;
  }


  /// Convenient function for streaming out lists of any streamable object.
  template<typename T>
  inline std::ostream& operator<<(std::ostream& os, const std::list<T>& vec) {
    os << "[ ";
    foreach (const T& i, vec) {
      os << i << " ";
    }
    os << "]";
    return os;
  }


}

#endif
