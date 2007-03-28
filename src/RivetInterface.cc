// -*- C++ -*-

#include "Rivet/Rivet.hh"
#include "Rivet/RivetHandler.hh"
#include "Rivet/Commandline.hh"
#include "Rivet/Tools/Logging.hh"
using namespace Rivet;

#include "HepMC/IO_Ascii.h"
#include "HepMC/GenEvent.h"
using namespace HepMC;


int main(int argc, char* argv[]) {
  Log& log = Log::getLog("Rivet.Main");

  // Configuration variables
  set<AnalysisName> cfgAnalyses;
  string cfgHistoFileName;
  HistoFormat cfgHistoFormat;
  map<string, Log::Level> cfgLogLevels;
  unsigned int cfgNumEvents;

  try {
      // Set up command line arguments
      TCLAP::CmdLine cmd("Analyse HepMC events from file using Rivet library", ' ', "0.1"); 
      //
      TCLAP::MultiArg<string>* analysesArg(0);
      TCLAP::SwitchArg* analysesAllArg(0);
      TCLAP::ValuesConstraint<string>* anaNameConstraint(0);
      Commandline::addAnalysisArgs(cmd, analysesArg, analysesAllArg, anaNameConstraint);
      //
      TCLAP::ValueArg<string>* histoNameArg(0);
      TCLAP::ValueArg<string>* histoTypeArg(0);
      TCLAP::ValuesConstraint<string>* histoTypeConstraint(0);
      Commandline::addHistoArgs(cmd, histoNameArg, histoTypeArg, histoTypeConstraint);
      //
      TCLAP::MultiArg<string>* logsArg(0);
      Commandline::addLoggingArgs(cmd, logsArg);
      //
      TCLAP::ValueArg<unsigned int> 
        numEventsArg("n", "numevents", "Max number of events to read (100000 by default)", 
                     false, 100000, "num", cmd);
      
      // Parse command line args
      cmd.parse(argc, argv);
      Commandline::useAnalysisArgs(cmd, analysesArg, analysesAllArg, cfgAnalyses);
      Commandline::useHistoArgs(cmd, histoNameArg, histoTypeArg, cfgHistoFileName, cfgHistoFormat);
      Commandline::useLoggingArgs(cmd, logsArg, cfgLogLevels);
      cfgNumEvents = numEventsArg.getValue();

  } catch (TCLAP::ArgException& e) { 
    cerr << "Command line error: " << e.error() << " for arg " << e.argId() << endl; 
    return EXIT_FAILURE;
  } catch (runtime_error& e) { 
    cerr << "Error: " << e.what() << endl; 
    return EXIT_FAILURE;
  }

  // Make a handler and add analyses
  RivetHandler rh;
  rh.addAnalyses(cfgAnalyses);
  rh.init();
  log << Log::INFO << "RivetHandler info: " << rh.info() << endl;

  // Make a HepMC output strategy
  /// @todo Add an option to specify an output file
  //HepMC::IO_Ascii hepmcOut("rivetgun.hepmc", std::ios::out);

  // Event loop
  HepMC::GenEvent myevent;
  for (unsigned int i = 0; i < cfgNumEvents; ++i) {    
    log << Log::INFO << "Event number " << i+1 << endl;
    rh.analyze(myevent);
  }
  log << Log::INFO << "Finished!"  << endl;

  // Finalise
  rh.finalize();

  return EXIT_SUCCESS;
}
