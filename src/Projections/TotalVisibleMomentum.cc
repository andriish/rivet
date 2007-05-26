// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/Cmp.hh"
#include "HepPDT/ParticleID.hh"

using namespace Rivet;


int TotalVisibleMomentum::compare(const Projection& p) const {
  const TotalVisibleMomentum & other = dynamic_cast<const TotalVisibleMomentum &>(p);
  return pcmp(*vfsproj, *other.vfsproj);
}


void TotalVisibleMomentum::project(const Event& e) {
  Log& log = getLog();

  _momentum = CLHEP::LorentzVector();
  _momentum.setPx(0.);
  _momentum.setPy(0.);
  _momentum.setPz(0.);
  _momentum.setE(0.);

  // Project into final state
  const VetoedFinalState& vfs = e.applyProjection(*vfsproj);

  // Get hadron and charge info for each particle, and fill counters appropriately
  for (ParticleVector::const_iterator p = vfs.particles().begin(); p != vfs.particles().end(); ++p) {
    //*_momentum+=p->getMomentum();
    _momentum+=p->getMomentum();
  }

  log << Log::DEBUG << "Done" << endl;

}


