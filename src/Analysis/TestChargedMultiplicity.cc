// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analysis/TestChargedMultiplicity.hh"
#include "HepPDT/ParticleID.hh"

using namespace Rivet;
using namespace HepMC;


TestChargedMultiplicity::~TestChargedMultiplicity() {}


void TestChargedMultiplicity::init() {
  /// @todo Book histogram here.
}


void TestChargedMultiplicity::analyze(const Event & event) {
  Logger& log = getLogger();
  log.setPriority(LogPriority::INFO);
  log << LogPriority::DEBUG << "Starting analyzing" << endlog;

  const Multiplicity& m = event.addProjection(mult);
  log << LogPriority::INFO << "Total multiplicity            = " << m.totalMultiplicity()           << endlog;
  log << LogPriority::INFO << "Total charged multiplicity    = " << m.totalChargedMultiplicity()    << endlog;
  log << LogPriority::INFO << "Total uncharged multiplicity  = " << m.totalUnchargedMultiplicity()  << endlog;
  log << LogPriority::INFO << "Hadron multiplicity           = " << m.hadronMultiplicity()          << endlog;
  log << LogPriority::INFO << "Hadron charged multiplicity   = " << m.hadronChargedMultiplicity()   << endlog;
  log << LogPriority::INFO << "Hadron uncharged multiplicity = " << m.hadronUnchargedMultiplicity() << endlog;

  /// @todo Fill histogram here

  log << LogPriority::DEBUG << "Finished analyzing" << endlog;
}


void TestChargedMultiplicity::finalize() {}


RivetInfo TestChargedMultiplicity::getInfo() const {
  return AnalysisBase::getInfo() + mult.getInfo();
}
