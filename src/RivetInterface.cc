// -*- C++ -*-

#include "Rivet/Projections/DISLepton.hh"

#include "Rivet/Rivet.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Tools/Commandline.hh"
#include "Rivet/Tools/Logging.hh"
using namespace Rivet;
using namespace TCLAP;

#include "HepMC/IO_Ascii.h"
#include "HepMC/GenEvent.h"
using namespace HepMC;


int main(int argc, char* argv[]) {

  // Configuration variables
  set<AnalysisName> cfgAnalyses;
  string cfgHistoFileName;
  HistoFormat cfgHistoFormat;
  Log::LevelMap cfgLogLevels;
  unsigned int cfgNumEvents;
  string cfgEventFile;

  try {
    // Set up command line arguments
    CmdLine cmd("Analyse HepMC events from a file using the Rivet system", ' ', "0.1"); 
    //
    MultiArg<string>* analysesArg(0);
    SwitchArg* analysesAllArg(0);
    ValuesConstraint<string>* anaNameConstraint(0);
    Commandline::addAnalysisArgs(cmd, analysesArg, analysesAllArg, anaNameConstraint);
    //
    ValueArg<string>* histoNameArg(0);
    ValueArg<string>* histoTypeArg(0);
    ValuesConstraint<string>* histoTypeConstraint(0);
    Commandline::addHistoArgs(cmd, histoNameArg, histoTypeArg, histoTypeConstraint);
    //
    MultiArg<string>* logsArg(0);
    Commandline::addLoggingArgs(cmd, logsArg);
    //
    ValueArg<unsigned int> 
      numEventsArg("n", "numevents", "Max number of events to read (100000 by default)", 
                   false, 100000, "num", cmd);
    
    UnlabeledValueArg<string>
      eventFileArg("eventfile", "File containing ASCII format HepMC events", true, "-", "filename", cmd);      
    
    // Parse command line args
    cmd.parse(argc, argv);
    Commandline::useAnalysisArgs(cmd, analysesArg, analysesAllArg, cfgAnalyses);
    Commandline::useHistoArgs(cmd, histoNameArg, histoTypeArg, cfgHistoFileName, cfgHistoFormat);
    Commandline::useLoggingArgs(cmd, logsArg, cfgLogLevels);
    cfgNumEvents = numEventsArg.getValue();
    cfgEventFile = eventFileArg.getValue();

  } catch (ArgException& e) { 
    cerr << "Command line error: " << e.error() << " for arg " << e.argId() << endl; 
    return EXIT_FAILURE;
  } catch (runtime_error& e) { 
    cerr << "Error: " << e.what() << endl; 
    return EXIT_FAILURE;
  }


  /////////////////////////////////////////////////////////////////////////////


  // Add all the log levels from the command line into the logging framework
  Log::setDefaultLevels(cfgLogLevels);
  Log& log = Log::getLog("Rivet.Main");


  // Make a handler and add analyses
  AnalysisHandler rh;
  for (AnalysisList::const_iterator ai = cfgAnalyses.begin(); ai != cfgAnalyses.end(); ++ai) {
     Analysis& a = Analysis::getAnalysis(*ai);
     BeamPair beams = a.getBeams();
     cout << a.name() << ": " << a.getBeams() << " " 
          << a.isCompatible(PROTON, PROTON) << endl;
     rh.addAnalysis(*ai);
  }
  rh.init();


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
    rh.analyze(myevent);
    ++num;
  }


  // Finish up
  log << Log::INFO << "Finished!"  << endl;
  rh.finalize();
  return EXIT_SUCCESS;
}
