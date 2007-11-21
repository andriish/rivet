// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Projections/CentralEtHCM.hh"
#include "Rivet/Cmp.hh"


namespace Rivet {

  int CentralEtHCM::compare(const Projection& p) const {
    const CentralEtHCM & other = dynamic_cast<const CentralEtHCM&>(p);
    return pcmp(_fshcm, other._fshcm);
  }


  void CentralEtHCM::project(const Event& e) {
    const FinalStateHCM& fs = e.applyProjection(_fshcm);
    _sumet = 0.0;
    for (int i=0, N=fs.particles().size(); i < N; ++i) {
      // Rapidity cut: abs rapidity < 0.5
      const FourMomentum p4 = fs.particles()[i].getMomentum();
      if (fabs(p4.rapidity()) < 0.5) _sumet += p4.Et();
    }
  }

}
