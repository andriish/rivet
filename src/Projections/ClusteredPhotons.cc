// -*- C++ -*-
#include "Rivet/Projections/ClusteredPhotons.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Cmp.hh"

namespace Rivet {

  int ClusteredPhotons::compare(const Projection& p) const {
    const PCmp fscmp = mkNamedPCmp(p, "FS");
    if (fscmp != PCmp::EQUIVALENT) return fscmp;

    const PCmp sigcmp = mkNamedPCmp(p, "Signal");
    if (sigcmp != PCmp::EQUIVALENT) return sigcmp;

    const ClusteredPhotons& other = dynamic_cast<const ClusteredPhotons&>(p);
    int rcmp = cmp(_dRmax, other._dRmax);
    return rcmp;
  }


  void ClusteredPhotons::project(const Event& e) {
    _theParticles.clear();
    if (!_dRmax>0.0) return;

    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    const FinalState& signal = applyProjection<FinalState>(e, "Signal");

    foreach (const Particle& p, fs.particles()) {
      bool clustered = false;
      if (p.pdgId() == PHOTON) {
        foreach (const Particle& l, signal.particles()) {
          FourMomentum p_l=l.momentum();
          FourMomentum p_P=p.momentum();
          if (deltaR(p_l.pseudorapidity(), p_l.azimuthalAngle(),
                     p_P.pseudorapidity(), p_P.azimuthalAngle()) < _dRmax) {
            clustered = true;
          }
        }
        if (clustered) _theParticles.push_back(p);
      }
    }
    getLog() << Log::DEBUG << name() <<" found " << _theParticles.size()
             <<" particles." << endl;
  }

}
