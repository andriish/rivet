// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Cmp.hh"
#include "Rivet/Tools/Utils.hh"
#include <algorithm>

namespace Rivet {

  int VetoedFinalState::compare(const Projection& p) const {
    const VetoedFinalState& other = dynamic_cast<const VetoedFinalState&>(p);
    const int fscmp = pcmp(_fsproj, other._fsproj);
    if (fscmp != 0) return fscmp;
    if (_vetoCodes == other._vetoCodes) return 0;
    if (_vetoCodes < other._vetoCodes) return -1; else return 1;
  }


  void VetoedFinalState::project(const Event& e) {
    Log log = getLog();
    const FinalState& fs = e.applyProjection(_fsproj);
    _theParticles.clear();
    _theParticles.reserve(fs.particles().size());
    const ParticleVector& fsps = fs.particles();
    for (ParticleVector::const_iterator p = fsps.begin(); p != fsps.end(); ++p) {
      if (log.isActive(Log::DEBUG)) {
        vector<long> codes; 
        for (VetoDetails::const_iterator code = _vetoCodes.begin(); code != _vetoCodes.end(); ++code) {
          codes.push_back(code->first);
        }
        const string codestr = "{ " + join(codes) + " }";
        log << Log::DEBUG << p->getPdgId() << " vs. veto codes = " 
            << codestr << " (" << codes.size() << ")" << endl;
      }
      const long pdgid = p->getPdgId();
      VetoDetails::iterator iter = _vetoCodes.find(pdgid);
      if ( (iter == _vetoCodes.end())) {
        log << Log::DEBUG << "Storing with PDG code " << pdgid << " pt " 
            << p->getMomentum().vector3().polarRadius() << endl;
        _theParticles.push_back(*p);
      } else {
        // This particle code is listed as a possible veto... check pT.
        BinaryCut ptrange = iter->second;
        // Make sure that the pT range is sensible.
        assert(ptrange.getHigherThan() <= ptrange.getLowerThan());
        double pt = p->getMomentum().vector3().polarRadius();
        stringstream rangess;
        if (ptrange.getHigherThan() < numeric_limits<double>::max()) rangess << ptrange.getHigherThan();
        rangess << " - ";
        if (ptrange.getLowerThan() < numeric_limits<double>::max()) rangess << ptrange.getLowerThan();
        log << Log::DEBUG << "ID = " << pdgid << ", pT range = " << rangess.str();
        stringstream debugline;
        debugline << "with PDG code = " << pdgid << " pT = " << p->getMomentum().vector3().polarRadius();
        if (pt < ptrange.getHigherThan() || pt > ptrange.getLowerThan()) { 
          log << Log::DEBUG << "Storing " << debugline << endl;
          _theParticles.push_back(*p);
        } else {
          log << Log::DEBUG << "Vetoing " << debugline << endl;
        }
      }
    }
  } 
  
}
