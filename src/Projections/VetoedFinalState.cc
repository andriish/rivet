// -*- C++ -*-

#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/Cmp.hh"
#include <algorithm>

using namespace std;

namespace Rivet {

  int VetoedFinalState::compare(const Projection& p) const {
    const VetoedFinalState& other = dynamic_cast<const VetoedFinalState&>(p);
    const int fscmp = pcmp(*_fsproj, *other._fsproj);
    if (fscmp != 0) return fscmp;
    if (_vetoCodes.size() != other._vetoCodes.size()) {
      return (_vetoCodes.size() < other._vetoCodes.size()) ? -1 : +1;
    } else {
      vector<long> myVetos = _vetoCodes;
      vector<long> otherVetos = other._vetoCodes;
      sort(myVetos.begin(), myVetos.end());
      sort(otherVetos.begin(), otherVetos.end());
      pair<vector<long>::iterator, vector<long>::iterator> wronguns =
        std::mismatch(myVetos.begin(), myVetos.end(), otherVetos.begin());
      if (wronguns.first == myVetos.end()) {
        return 0;
      } else {
        return (*(wronguns.first) < *(wronguns.second)) ? -1 : +1;
      }
    }
  }

  void VetoedFinalState::project(const Event& e) {
    Log log = getLog();
    const FinalState& fs = e.applyProjection(*_fsproj);
    _theParticles.clear();
    _theParticles.reserve(fs.particles().size());
    const ParticleVector fsps = fs.particles();
    for (ParticleVector::const_iterator p = fsps.begin(); p != fsps.end(); ++p) {
      log << Log::DEBUG << p->getPdgId() << " codes " << _vetoCodes.size() << endl;
      log << Log::DEBUG <<_vetoCodes.at(0) << "," << _vetoCodes.at(1) << "," <<_vetoCodes.at(2) << "," <<_vetoCodes.at(3) << "," <<_vetoCodes.at(4) << "," <<_vetoCodes.at(5) << "," <<_vetoCodes.at(6) << "," << endl;
      if (find(_vetoCodes.begin(), _vetoCodes.end(), p->getPdgId()) == _vetoCodes.end()) {
	log << Log::DEBUG << "stored" << endl;
        _theParticles.push_back(*p);
      }
    }
  }

}
