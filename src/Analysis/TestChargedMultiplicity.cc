// -*- C++ -*-

#include "Rivet/Analysis/TestChargedMultiplicity.h"
#include "Rivet/Logging.h"
#include "HepPDT/ParticleID.hh"

using namespace Rivet;
using namespace std;
using namespace HepMC;


TestChargedMultiplicity::~TestChargedMultiplicity() {}


void TestChargedMultiplicity::init() {
  /// @todo Book histogram here.
}


void TestChargedMultiplicity::analyze(const Event & event) {
  Logger& log = getLogger();
  log.setPriority(LogPriority::INFO);
  log << LogPriority::DEBUG << "Starting analyzing" << endlog;

  const FinalStateProjection& fs = event.addProjection(fsproj);
  unsigned int chmult(0), unchmult(0);
  unsigned int particleNum(0);
  for (PVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
    ++particleNum;
    HepPDT::ParticleID pInfo(p->id);
    if (pInfo.isHadron()) {
      if (pInfo.threeCharge() != 0) {
        ++chmult;
        log << LogPriority::DEBUG << "Incrementing charged multiplicity   = " << chmult << endlog;
      } else {
        ++unchmult;
        log << LogPriority::DEBUG << "Incrementing uncharged multiplicity = " << unchmult << endlog;
      }
    }
  }
  log << LogPriority::INFO << "Event charged multiplicity   = " << chmult << endlog;
  log << LogPriority::INFO << "Event uncharged multiplicity = " << unchmult << endlog;

  /// @todo Fill histogram here

  log << LogPriority::DEBUG << "Finished analyzing" << endlog;
}


void TestChargedMultiplicity::finalize() {}


RivetInfo TestChargedMultiplicity::getInfo() const {
  return AnalysisBase::getInfo() + fsproj.getInfo();
}
