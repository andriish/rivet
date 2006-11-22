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



  /// @todo Fill histogram here

  log << LogPriority::DEBUG << "Finished analyzing" << endlog;
}


void TestChargedMultiplicity::finalize() {}


RivetInfo TestChargedMultiplicity::getInfo() const {
  return AnalysisBase::getInfo() + fsproj.getInfo();
}
