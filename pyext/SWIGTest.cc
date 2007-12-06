// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/AnalysisLoader.hh"
//#include "Rivet/Tools/Logging.hh"
using namespace Rivet;

#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
using namespace HepMC;

using namespace std;


int run() {

  // Make a handler and add analyses
  //AnalysisHandler handler;
  //for (set<string>::const_iterator an = cfgAnalyses.begin(); an != cfgAnalyses.end(); ++an) {
  Analysis* a = AnalysisLoader::getAnalysis("EXAMPLE");
  //BeamPair beams = a->getBeams();
  cout << "Analysis name: " << a->getName() << endl;
  //log << Log::INFO << "Cuts:" << a->getCuts() << endl;
  //handler.addAnalysis(a->getName());
  delete a;
  //}
  //handler.init();

  // Finish up
  cout << "Finished!"  << endl;
  //handler.finalize();
  return EXIT_SUCCESS;
}
