#ifndef RIVET_RIVETBOOST_HH
#define RIVET_RIVETBOOST_HH

#include "boost/smart_ptr.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/assign.hpp"

#include "boost/foreach.hpp"
#define foreach BOOST_FOREACH

namespace Rivet {
  using boost::shared_ptr;
  using boost::lexical_cast;
}

#endif
