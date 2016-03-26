#ifndef RIVET_RIVETBOOST_HH
#define RIVET_RIVETBOOST_HH

/// @todo Replaced by C++11 std?
// #include "boost/smart_ptr.hpp"
// #include "boost/pointer_cast.hpp"
// #include "boost/type_traits.hpp"
// #include "boost/utility.hpp"

/// @todo Not in C++11 -- can be removed / mocked up easily?
// #include "boost/lexical_cast.hpp"
// #include "boost/algorithm/string.hpp"
#include "boost/assign.hpp"

namespace Rivet {


  // // Smart pointers
  // using boost::shared_ptr;

  // Clever casts
  /// @todo Replace?
  // using boost::lexical_cast;
  // using boost::bad_lexical_cast;

  // Clever assignment shortcuts
  /// @todo Replace with universal brace init
  using namespace boost::assign;

  // General Boost namespace
  /// @todo Remove?
  // using namespace boost;


}

#endif
