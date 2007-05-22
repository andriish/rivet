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
      //cout << LESS_EQ << endl;
      if (value < _cuts[quantity].lowerthan()) {
        _cuts[quantity].lowerthan() = value;
      }
      break;
    case MORE_EQ:
      //cout << MORE_EQ << endl;
      if (value > _cuts[quantity].higherthan()) {
        _cuts[quantity].higherthan() = value;
      }
      break;
    case EQUAL:
      //cout << EQUAL << endl;
      _cuts[quantity].lowerthan() = value;
      _cuts[quantity].higherthan() = value;
      break;
    }
    //cout << "Cuts::addCut(): " << quantity << " <= " << _cuts[quantity].lowerthan() << " (" << value << ")" << endl;
    //cout << "Cuts::addCut(): " << quantity << " >= " << _cuts[quantity].higherthan() << " (" << value << ")" << endl;

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
      os << endl;
      os << setw(12) << std::left << cut->first;
      if (cut->second.higherthan() > -numeric_limits<double>::max()) {
        os << setw(3) << ">=";
        os << setw(10) << std::right << cut->second.higherthan();
      } else {
        os << setw(13) << "";
      }
      if (cut->second.lowerthan() < numeric_limits<double>::max()) {
        os << setw(3) << "<=";
        os << setw(10) << std::right << cut->second.lowerthan();
      } else {
        os << setw(13) << "";
      }
    }
    return os;
  }


}
