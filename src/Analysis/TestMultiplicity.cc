// -*- C++ -*-

#include "Rivet/Analysis/TestMultiplicity.hh"

using namespace Rivet;
using namespace std;


TestMultiplicity::~TestMultiplicity() {}


void TestMultiplicity::init() {
  // Book histogram here.
}


void TestMultiplicity::analyze(const Event & event) {
  Logger& log = getLogger();
  log.setPriority(LogPriority::INFO);

  const FinalStateProjection& fs = event.addProjection(fsproj);
  int mult =  fs.particles().size();
  log << LogPriority::INFO << "Event multiplicity = " << mult << endlog;

  // Fill histogram here.
}


void TestMultiplicity::finalize() {}


RivetInfo TestMultiplicity::getInfo() const {
  return AnalysisBase::getInfo() + fsproj.getInfo();
}
