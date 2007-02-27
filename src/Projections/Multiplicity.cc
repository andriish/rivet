// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Multiplicity.hh"
#include "Rivet/Projections/Cmp.hh"
#include "HepPDT/ParticleID.hh"

using namespace Rivet;


int Multiplicity::compare(const Projection & p) const {
  const Multiplicity & other = dynamic_cast<const Multiplicity &>(p);
  return pcmp(*fsproj, *other.fsproj);
}


void Multiplicity::project(const Event & e) {
  Logger& log = getLogger();
  log.setPriority(LogPriority::INFO);

  // Clear counters
  totalMult_ = 0;
  totalChMult_ = 0;
  totalUnchMult_ = 0;
  hadMult_ = 0;
  hadChMult_ = 0;
  hadUnchMult_ = 0;

  // Project into final state
  const FinalState& fs = e.applyProjection(*fsproj);

  // Get hadron and charge info for each particle, and fill counters appropriately
  for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
    ++totalMult_;
    HepPDT::ParticleID pInfo = p->id;
    bool isHadron = pInfo.isHadron();
    if (pInfo.threeCharge() != 0) {
      ++totalChMult_;
      if (isHadron) {
        ++hadMult_;
        ++hadChMult_;
      }
      log << LogPriority::DEBUG << "Incrementing charged multiplicity   = " << totalChMult_ 
          << " (" << hadChMult_ << " hadrons)" << endlog;
    } else {
      ++totalUnchMult_;
      if (isHadron) {
        ++hadMult_;
        ++hadUnchMult_;
      }
      log << LogPriority::DEBUG << "Incrementing uncharged multiplicity = " << totalUnchMult_
          << " (" << hadUnchMult_ << " hadrons)" << endlog;
    }
  }
}

RivetInfo Multiplicity::getInfo() const {
  return Projection::getInfo() + fsproj->getInfo();
}

