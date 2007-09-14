// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Multiplicity.hh"
#include "Rivet/Cmp.hh"
#include "HepPDT/ParticleID.hh"


namespace Rivet {

  int Multiplicity::compare(const Projection& p) const {
    const Multiplicity& other = dynamic_cast<const Multiplicity&>(p);
    return pcmp(*_fsproj, *other._fsproj);
  }


  void Multiplicity::project(const Event& e) {
    Log& log = getLog();

    // Clear counters
    _totalMult = 0;
    _totalChMult = 0;
    _totalUnchMult = 0;
    _hadMult = 0;
    _hadChMult = 0;
    _hadUnchMult = 0;

    // Project into final state
    const FinalState& fs = e.applyProjection(*_fsproj);

    // Get hadron and charge info for each particle, and fill counters appropriately
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      ++_totalMult;
      HepPDT::ParticleID pInfo = p->getPdgId();
      bool isHadron = pInfo.isHadron();
      if (pInfo.threeCharge() != 0) {
        ++_totalChMult;
        if (isHadron) {
          ++_hadMult;
          ++_hadChMult;
        }
        log << Log::DEBUG << "Incrementing charged multiplicity = " << _totalChMult
            << " (" << _hadChMult << " hadrons)" << endl;
      } else {
        ++_totalUnchMult;
        if (isHadron) {
          ++_hadMult;
          ++_hadUnchMult;
        }
        log << Log::DEBUG << "Incrementing uncharged multiplicity = " << _totalUnchMult
            << " (" << _hadUnchMult << " hadrons)" << endl;
      }
    }
  }


}
