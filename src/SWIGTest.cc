// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Tools/Commandline.hh"
#include "Rivet/Tools/Logging.hh"
using namespace Rivet;

#include "HepMC/IO_Ascii.h"
#include "HepMC/GenEvent.h"
using namespace HepMC;


int runRivet() {

  // Configuration variables
  set<string> cfgAnalyses;
  cfgAnalyses.insert("TEST");
  size_t cfgNumEvents = 10;
  string cfgEventFile = "test.hepmc";
  Log& log = Log::getLog("Rivet.Main");


  // Make a handler and add analyses
  AnalysisHandler handler;
  for (set<string>::const_iterator an = cfgAnalyses.begin(); an != cfgAnalyses.end(); ++an) {
    Analysis* a = AnalysisLoader::getAnalysis(*an);
    BeamPair beams = a->getBeams();
    log << Log::INFO << "Analysis name: " << a->getName() << endl;
    log << Log::DEBUG << a->getBeams() << a->isCompatible(PROTON, PROTON) << endl;
    log << Log::INFO << "Cuts:" << a->getCuts() << endl;
    handler.addAnalysis(a->getName());
    delete a;
  }
  handler.init();


  // Make a HepMC input
  IO_Ascii hepmcIn(cfgEventFile.c_str(), std::ios::in);
  if (hepmcIn.rdstate() != 0) {
    log << Log::ERROR << "Couldn't read HepMC events from file: " << cfgEventFile << endl;
  }


  // Event loop
  unsigned int num = 0;
  while (num < cfgNumEvents) {
    HepMC::GenEvent myevent;
    if (! hepmcIn.fill_next_event(&myevent)) break;
    log << Log::INFO << "Event number " << num + 1 << endl;
    handler.analyze(myevent);
    ++num;
  }


  // Finish up
  log << Log::INFO << "Finished!"  << endl;
  handler.finalize();
  return 0;
}
