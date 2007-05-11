// -*- C++ -*-

#include "Rivet/Cut.hh"

using namespace Rivet;
using namespace std;


namespace Rivet {


  Cuts& Cuts::addCut(const string& quantity, const Comparison& comparison, const double value) {
    // If this quantity doesn't yet have any associated cuts, make a
    // default cut with effectively infinitely loose cuts.
    if (_cuts.find(quantity) == _cuts.end()) {
      _cuts[quantity] = BinaryCut();
    }
    // Combine cuts in the most restrictive way.
    switch (comparison) {
    case LESS_EQ: 
      if (value < _cuts[quantity].lowerthan()) 
        _cuts[quantity].lowerthan() = value;
    case GREATER_EQ: 
      if (value > _cuts[quantity].higherthan()) 
        _cuts[quantity].higherthan() = value;
    case EQUAL: 
      _cuts[quantity].lowerthan() = value;
      _cuts[quantity].higherthan() = value;
    }
    // Allow method chaining.
    return *this;
  }



  bool Cuts::checkConsistency() const { 
    for (Cuts::const_iterator c = begin(); c != end(); ++c) {
      if (c->second.lowerthan() < c->second.lowerthan()) {
        ostringstream msg;
        msg << "Constraints on " << c->first << " are incompatible: "
            << ">=" << c->second.higherthan() << " AND "
            << "<=" << c->second.lowerthan();
        throw runtime_error(msg.str());
      }
    }
    return true;
  }



  ostream& Cuts::print(ostream & os) const {
    for (Cuts::const_iterator cut = begin(); cut != end(); ++cut) {
      os << std::left;
      os << setw(12) << cut->first;
      if (cut->second.higherthan() > numeric_limits<double>::min()) {
        os << setw(10) << ">=" << cut->second.higherthan();
      } else {
        os << setw(10) << "";
      }
      if (cut->second.lowerthan() < numeric_limits<double>::max()) {
        os << setw(10) << "<=" << cut->second.lowerthan();
      } else {
        os << setw(10) << "";
      }
      os << endl;
    }
    return os;
  }


}
