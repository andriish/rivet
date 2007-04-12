// -*- C++ -*-

#include "Rivet/RivetInfo.hh"
#include <sstream>
#include <cmath>

using namespace Rivet;
using namespace std;


RivetInfo& RivetInfo::addParam(Quantity quantity, Comparison comparison, double value) {
  ParamKey key(quantity, comparison);
  Params::iterator p = _params.find(key);
  if (p == _params.end()) {
    _params[key] = value;
  } else {
    switch (comparison) {
    case LESS: 
      if (value < p->second) _params[key] = value;
    case GREATER: 
      if (value > p->second) _params[key] = value;
    case EQUAL:
      _params[key] = value;
    }
  }
  return *this;
}


RivetInfo& RivetInfo::addValidBeamPair(BeamParticle b1, BeamParticle b2, double mom1, double mom2) {
  //cout << _validBeams.size() << endl;
  Beam beam1(b1, mom1);
  Beam beam2(b2, mom2);
  BeamPair beampair(beam1, beam2);
  _validBeams.insert(beampair);
  //cout << _validBeams.size() << endl;
  _resolveValidBeams();
  return *this;
}


void RivetInfo::_check() const { 
  for (Params::const_iterator p1 = _params.begin(); p1 != _params.end(); ++p1) {
    ParamKey key1 = p1->first;
    for (Params::const_iterator p2 = p1; p2 != _params.end(); ++p2) {
      ParamKey key2 = p2->first;
      if (key1.first == key2.first) {
        if ((key1.second == LESS && key2.second == GREATER && p1->second <= p2->second) ||
            (key2.second == LESS && key1.second == GREATER && p2->second <= p1->second)) {
          ostringstream msg;
          msg << "Constraints on " << quantityToStr(key1.first) << " are incompatible: "
              << comparisonToStr(key1.second) << p1->second << " AND "
              << comparisonToStr(key2.second) << p2->second;
          throw runtime_error(msg.str());
        }
      }
    }
  }
}


void RivetInfo::_combineParams(const Params& otherps) {
  for (Params::const_iterator op = otherps.begin(); op != otherps.end(); ++op) {
    addParam(op->first.first, op->first.second, op->second);
  }
}


bool same(BeamParticle a, BeamParticle b) {
  return (a == b || a == ANY || b == ANY);
}

void RivetInfo::_combineValidBeams(const BeamsSet& otherbs) {
  for (BeamsSet::const_iterator bs = otherbs.begin(); bs != otherbs.end(); ++bs) {
    _validBeams.insert(*bs);
  }
  _resolveValidBeams();
}

void RivetInfo::_resolveValidBeams() {
  BeamsSet myBeams(_validBeams);
  _validBeams.clear();
  BeamsSet::const_iterator bs1, bs2;

  for (BeamsSet::const_iterator bs1 = myBeams.begin(); bs1 != myBeams.end(); ++bs1) {
    for (BeamsSet::const_iterator bs2 = otherbs.begin(); bs2 != otherbs.end(); ++bs2) {
      // Make pairs of candidates to be particle 1 and particle 2
      Beam b1a(bs1->first), b2a(bs1->second), b1b, b2b;

      // If the particle names match...
      if (same(bs1->first.first, bs2->first.first) && same(bs1->second.first, bs2->second.first)) {
        b1b = bs2->first;
        b2b = bs2->second;
      } else if (same(bs1->first.first, bs2->second.first) && same(bs1->second.first, bs2->first.first)) {
        b1b = bs2->second;
        b2b = bs2->first;
      } else {
        // These beams don't match, so don't add them to the new set
        continue;
      }
      // Choose the more restrictive member of each pairing
      BeamParticle p1 = (b1a.first == ANY) ? b1b.first : b1a.first;
      BeamParticle p2 = (b2a.first == ANY) ? b2b.first : b2a.first;

      // Now deal with the momenta
      double mom1 = -1.0;
      if (b1a.second >= 0 && b1b.second >= 0 && fabs(b1a.second - b1b.second) > 1e-3) 
        throw runtime_error("Matching valid beams have incompatible allowed momenta");
      else mom1 = (b1a.second < b1b.second) ? b1a.second : b1b.second;
      double mom2 = -1.0;
      if (b2a.second >= 0 && b2b.second >= 0 && fabs(b2a.second - b2b.second) > 1e-3) 
        throw runtime_error("Matching valid beams have incompatible allowed momenta");
      else mom2 = (b2a.second < b2b.second) ? b2a.second : b2b.second;

      // Add combined beams
      addValidBeamPair(p1, p2, mom1, mom2);
    }
  }
}


ostream& RivetInfo::print(ostream & os) const {
  os << "Rivet info:" << endl;
  os << "# params = " << _params.size() << endl;
  os << std::left;
  for (Params::const_iterator p = _params.begin(); p != _params.end(); ++p) {
    os << setw(12) << quantityToStr(p->first.first)
       << setw(12) << comparisonToStr(p->first.second)
       << setw(10) << p->second 
       << endl;
  }
  os << "# beam pairs = " << _validBeams.size() << endl;
  for (BeamsSet::const_iterator b = _validBeams.begin(); b != _validBeams.end(); ++b) {
    os << "[" 
       << getKnownBeamParticles()[b->first.first] << "(" << b->first.second << ")"
       << getKnownBeamParticles()[b->second.first] << "(" << b->second.second << ")"
       << "]";
  }
  return os;
}

