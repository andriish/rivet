// -*- C++ -*-

#include "Rivet/Cut.hh"
#include <sstream>
#include <cmath>

using namespace Rivet;
using namespace std;


// RivetInfo& RivetInfo::addParam(Quantity quantity, Comparison comparison, double value) {
//   ParamKey key(quantity, comparison);
//   Params::iterator p = _params.find(key);
//   if (p == _params.end()) {
//     _params[key] = value;
//   } else {
//     switch (comparison) {
//     case LESS: 
//       if (value < p->second) _params[key] = value;
//     case GREATER: 
//       if (value > p->second) _params[key] = value;
//     case EQUAL:
//       _params[key] = value;
//     }
//   }
//   return *this;
// }


// void RivetInfo::_check() const { 
//   for (Params::const_iterator p1 = _params.begin(); p1 != _params.end(); ++p1) {
//     ParamKey key1 = p1->first;
//     for (Params::const_iterator p2 = p1; p2 != _params.end(); ++p2) {
//       ParamKey key2 = p2->first;
//       if (key1.first == key2.first) {
//         if ((key1.second == LESS && key2.second == GREATER && p1->second <= p2->second) ||
//             (key2.second == LESS && key1.second == GREATER && p2->second <= p1->second)) {
//           ostringstream msg;
//           msg << "Constraints on " << quantityToStr(key1.first) << " are incompatible: "
//               << comparisonToStr(key1.second) << p1->second << " AND "
//               << comparisonToStr(key2.second) << p2->second;
//           throw runtime_error(msg.str());
//         }
//       }
//     }
//   }
// }


// void RivetInfo::_combineParams(const Params& otherps) {
//   for (Params::const_iterator op = otherps.begin(); op != otherps.end(); ++op) {
//     addParam(op->first.first, op->first.second, op->second);
//   }
// }


// bool same(BeamParticle a, BeamParticle b) {
//   return (a == b || a == ANY || b == ANY);
// }


ostream& Cut::print(ostream & os) const {
  os << std::left;
  os << setw(12) << toString(_quantity)
     << setw(12) << toString(_comp)
     << setw(10) << _value 
     << endl;
  return os;
}

