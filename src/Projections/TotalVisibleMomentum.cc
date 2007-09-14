// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Cmp.hh"
#include "HepPDT/ParticleID.hh"


namespace Rivet {

  int TotalVisibleMomentum::compare(const Projection& p) const {
    const TotalVisibleMomentum& other = dynamic_cast<const TotalVisibleMomentum&>(p);
    return pcmp(*_fsproj, *other._fsproj);
  }


  void TotalVisibleMomentum::project(const Event& e) {
    Log& log = getLog();

    _momentum = LorentzVector();
    _momentum.setPx(0.0);
    _momentum.setPy(0.0);
    _momentum.setPz(0.0);
    _momentum.setE(0.0);

    // Project into final state
    const FinalState& fs = e.applyProjection(*_fsproj);

    // Get hadron and charge info for each particle, and fill counters appropriately
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      _momentum += p->getMomentum();
    }

    log << Log::DEBUG << "Done" << endl;
  }


}
