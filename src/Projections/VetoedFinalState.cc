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
    if (_vetoCodes == other._vetoCodes) return 0;
    if (_vetoCodes < other._vetoCodes) return -1; else return 1;
  }


  void VetoedFinalState::project(const Event& e) {
    Log log = getLog();
    const FinalState& fs = e.applyProjection(*_fsproj);
    _theParticles.clear();
    _theParticles.reserve(fs.particles().size());
    const ParticleVector& fsps = fs.particles();
    for (ParticleVector::const_iterator p = fsps.begin(); p != fsps.end(); ++p) {
      if (log.isActive(Log::DEBUG)) {
        stringstream codes; 
        codes << "{";
        for (VetoDetails::const_iterator code = _vetoCodes.begin(); code != _vetoCodes.end(); ++code) {
          codes << " " << code->first;
        }
        codes << " }";
        log << Log::DEBUG << p->getPdgId() << " vs. veto codes = " 
            << codes.str() << " (" << _vetoCodes.size() << ")" << endl;
      }
      const long pdgid = p->getPdgId();
      VetoDetails::iterator iter = _vetoCodes.find(pdgid);
      if ( (iter == _vetoCodes.end())) {
        log << Log::DEBUG << "Storing with PDG code " << pdgid << " pt " << p->getMomentum().perp() << endl;
        _theParticles.push_back(*p);
      } else {
        // This particle code is listed as a possible veto... check pT.
        pair<double, double> ptrange = iter->second;
        // Make sure that the pT range is sensible.
        assert(ptrange.first <= ptrange.second);
        double pt = p->getMomentum().perp();
        stringstream rangess;
        if (ptrange.first < numeric_limits<double>::max()) rangess << ptrange.first;
        rangess << " - ";
        if (ptrange.second < numeric_limits<double>::max()) rangess << ptrange.second;
        log << Log::DEBUG << "ID = " << pdgid << ", pT range = " << rangess.str();
        stringstream debugline;
        debugline << "with PDG code = " << pdgid << " pT = " << p->getMomentum().perp();
        if (pt < ptrange.first || pt > ptrange.second) { /// @todo Is this the right way round?
          log << Log::DEBUG << "Storing " << debugline << endl;
          _theParticles.push_back(*p);
        } else {
          log << Log::DEBUG << "Vetoing " << debugline << endl;
        }
      }
    }
  } 
  
}
