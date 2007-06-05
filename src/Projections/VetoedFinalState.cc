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
        for (set<long>::const_iterator code = _vetoCodes.begin(); 
             code != _vetoCodes.end(); ++code) codes << " " << *code;
        codes << " }";
        log << Log::DEBUG << p->getPdgId() << " vs. veto codes = " 
            << codes.str() << " (" << _vetoCodes.size() << ")" << endl;
      }
      const long pdgid = p->getPdgId();
      if (_vetoCodes.find(pdgid) == _vetoCodes.end()) {
        log << Log::DEBUG << "Storing with PDG code " << pdgid << endl;
        _theParticles.push_back(*p);
        //cout << "VetoedFinalState.cc: Particle (pdgid=" << pdgid << " accepted, _vetoCodes.size()=" << 
        //_vetoCodes.size() << "   eta=" << p->getMomentum().eta()  << endl;
      }
      //else {cout << "VetoedFinalState.cc: p->getPdgId()=" << pdgid << ": particle skipped!" << endl; }
    }
  }


}
