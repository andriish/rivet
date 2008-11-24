// -*- C++ -*-

#include "Rivet/Config/BuildOptions.hh"

#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/Configuration.hh"
using namespace Rivet;

#ifdef HAVE_AGILE
#include "AGILe/AGILe.hh"
#include "AGILe/Particle.hh"
#include "AGILe/Generator.hh"
#include "AGILe/Loader.hh"
using namespace AGILe;
#endif

#include "HepMC/IO_GenEvent.h"
using namespace HepMC;

#include <csignal>


namespace {

  volatile sig_atomic_t RECVD_KILL_SIGNAL = 0;

  extern "C" {
    void handleKillSignal(int sigid) {
      RECVD_KILL_SIGNAL = sigid;
      signal(sigid, SIG_DFL);
    }
  }

  void registerKillSignalHandler() {
    signal(SIGTERM, handleKillSignal);
    signal(SIGHUP,  handleKillSignal);
    signal(SIGINT,  handleKillSignal);
    signal(SIGUSR2, handleKillSignal);
  }

}



namespace Rivet {


  #ifdef HAVE_AGILE
  void setupGenerator(Configuration& cfg, Log& log, Generator*& gen) {
    // Load generator libraries
    Loader::initialize();
    log << Log::INFO << "Requested generator = " << cfg.generatorName << endl;
    try {
      Loader::loadGenLibs(cfg.generatorName);
    } catch (Error& e) {
      log << Log::ERROR << "Error when loading the generator library:"
          << endl << e.what() << endl;
      throw Error("Error when loading the generator library"); 
    }
    
    gen = Loader::createGen();
    if (!gen) {
      log << Log::ERROR << "No generator chosen... exiting" << endl;
      throw Error("No generator chosen");
    }
    // Seed random number generator
    gen->setSeed(cfg.rngSeed);
    log << Log::DEBUG << "Finished setting seed " << endl;
    
    // Param passing
    for (ParamMap::const_iterator p = cfg.params.begin(); p != cfg.params.end(); ++p) {
      log << Log::INFO << "Setting param " << p->first << ": " << p->second << endl;
      gen->setParam(p->first, p->second);
    }
    log << Log::DEBUG << "Finished setting parameters " << endl;
    
    // Set initial state
    try {
      gen->setInitialState(cfg.beam1, cfg.mom1, cfg.beam2, cfg.mom2);
    } catch (Error& e) {
      log << Log::ERROR << "Beam particle error: " << cfg.beam1 << " or " 
          << cfg.beam2 << " is invalid" << endl;
      throw Error("Invalid beam particle");
    }
    log << Log::DEBUG << "Finished setting initial state " << endl;
  } 
  #endif


  void setupRivet(Configuration& cfg, Log& log, AnalysisHandler& rh) {
    // Initialise Rivet
    if (cfg.runRivet) {
      vector<string> analyses;
      foreach (const string& aname, cfg.analyses) {
        Analysis* a = AnalysisLoader::getAnalysis(aname);
        if (a == 0) {
          log << Log::WARN << "Cannot find an analysis to match name `" << aname << "'" << endl;
          log << Log::WARN << "If `" << aname << "' is a user analysis, is it " 
              << "registered in the `getAnalysisBuilders' function as`" << toUpper(aname) << "?" << endl;
          break;
        }
        const bool beamMatch = a->isCompatible(cfg.beam1, cfg.beam2);
        log << Log::DEBUG << a->name() << ": " << a->requiredBeams() << " "
            << (beamMatch ? "MATCH" : "NO-MATCH") << endl;
        if (beamMatch || cfg.readHepMC) {
          log << Log::INFO << "Using analysis " << a->name() << endl;
          analyses.push_back(aname);
        } else {
          log << Log::WARN << "Removing inappropriate analysis " << a->name() << endl;
          delete a;
        }
      }
      /// @todo Decide who should be controlling the proj handler (do it via AnalysisHandler?)
      ProjectionHandler::create()->clear();
      rh.addAnalyses(analyses);
      rh.init();
    }
  }


  // Notify about event number
  void printEventNumber(long num) {
    Log& nevtlog = Log::getLog("RivetGun.NEvt");
    Log::Level lev = Log::DEBUG;
    if (round(num/100.0) == num/100.0) lev = Log::INFO;
    if (round(num/1000.0) == num/1000.0) lev = Log::WARN;
    nevtlog << lev << "Event number " << num << endl;
  }


  void useEvent(Configuration& cfg, const HepMC::GenEvent& myevent, 
                AnalysisHandler& rh, IO_GenEvent* hepmcOut) {
    // Print out HepMC event if in DEBUG mode
    Log& hepmclog = Log::getLog("RivetGun.HepMC");
    if (hepmclog.isActive(Log::DEBUG)) {
      myevent.print();
    }
      
    // Run Rivet analyses
    if (cfg.runRivet) rh.analyze(myevent);
    
    // Write out event to file
    if (cfg.writeHepMC) {
      hepmcOut->write_event(&myevent);
      //myevent.print(cout);
    }
  }



  #ifdef HAVE_AGILE
  void doGenEventLoop(Configuration& cfg, Log& log, Generator* gen, 
                      AnalysisHandler& rh, IO_GenEvent* hepmcOut) {
    log << Log::INFO << "Generating " << cfg.numEvents << " events." << endl;
    HepMC::GenEvent myevent;

    //Rivet::Beam beamproj;

    for (size_t i = 0; i < cfg.numEvents; ++i) {    
      assert(gen);
      gen->makeEvent(myevent);
      printEventNumber(i+1);
      useEvent(cfg, myevent, rh, hepmcOut);

      // Exit nicely if we've been signalled during this iteration
      if (RECVD_KILL_SIGNAL != 0) {
        log << Log::WARN << "Leaving event loop early due to signal " 
            << RECVD_KILL_SIGNAL << endl;
        break;
      }
    }
    log << Log::INFO << "Finished generator event loop!"  << endl;
  }
  #endif



  void doReadEventLoop(Configuration& cfg, Log& log, IO_GenEvent* hepmcIn, 
                       AnalysisHandler& rh, IO_GenEvent* hepmcOut) {
    log << Log::INFO << "Reading " << cfg.numEvents << " events." << endl;
    HepMC::GenEvent myevent;
    for (size_t i = 0; i < cfg.numEvents; ++i) {
      // Check stream state
      const int readstate = hepmcIn->rdstate();
      if (readstate != 0) {
        log << Log::ERROR << "Bad HepMC input stream state: " << cfg.hepmcInFile
            << " -> " << readstate << endl;
        break;
      }
      
      // Fill and use next event
      myevent.clear();
      const bool gotevent = hepmcIn->fill_next_event(&myevent);
      if (gotevent) {
        printEventNumber(i+1);
        useEvent(cfg, myevent, rh, hepmcOut);
      } else {
        log << Log::ERROR << "Couldn't read next HepMC event from file: " << cfg.hepmcInFile << endl;
        break;
      }

      // Exit nicely if we've been signalled during this iteration
      if (RECVD_KILL_SIGNAL != 0) {
        log << Log::WARN << "Leaving event loop early due to signal " 
            << RECVD_KILL_SIGNAL << endl;
        break;
      }
    }
    log << Log::INFO << "Finished read-in event loop!"  << endl;
  }



  void generate(Configuration& cfg, Log& log) {
    log << Log::DEBUG << "Entered generate() function " << endl;

    // Make sure that the event loop can exit gracefully
    registerKillSignalHandler();

    #ifdef HAVE_AGILE 
    Generator* gen = 0;
    if (!cfg.readHepMC) { 
      setupGenerator(cfg, log, gen);
    }
    #endif

    // Configure analyses etc
    AnalysisHandler rh(cfg.histoName, cfg.histoFormat);
    setupRivet(cfg, log, rh);

    // Make a HepMC output
    IO_GenEvent* hepmcOut(0);
    if (cfg.writeHepMC) {
      log << Log::DEBUG << "Making a HepMC output stream to " << cfg.hepmcOutFile << endl;
      hepmcOut = new IO_GenEvent(cfg.hepmcOutFile.c_str(), std::ios::out);
    }

    // Make a HepMC input
    IO_GenEvent* hepmcIn(0);
    if (cfg.readHepMC) {
      log << Log::DEBUG << "Making a HepMC input stream from " << cfg.hepmcInFile << endl;
      hepmcIn = new IO_GenEvent(cfg.hepmcInFile.c_str(), std::ios::in);
      log << Log::DEBUG << "Testing HepMC input stream from " << cfg.hepmcInFile << endl;
      if (hepmcIn && hepmcIn->rdstate() != 0) {
        log << Log::ERROR << "Couldn't read HepMC event file: " << cfg.hepmcInFile << endl;
        throw Error("Couldn't read HepMC event file: " + cfg.hepmcInFile);
      }
    }
    
    // Run event loop
    if (cfg.readHepMC) {
      doReadEventLoop(cfg, log, hepmcIn, rh, hepmcOut);
    } else {
      #ifdef HAVE_AGILE
      doGenEventLoop(cfg, log, gen, rh, hepmcOut);
      #else
      throw Error("RivetGun not compiled with AGILe: events can only be read from HepMC files");
      #endif
    }

    // Finalise Rivet
    #ifdef HAVE_AGILE
    if (gen) gen->finalize();
    #endif
    delete hepmcOut;
    delete hepmcIn;
    if (cfg.runRivet) {
      if (rh.needCrossSection()) {
        #ifdef HAVE_AGILE
        if (gen) {
          rh.setCrossSection(gen->getCrossSection());
        } else {
          log << Log::ERROR << "Warning: Cross section needed but no Generator created." <<endl;
          log << Log::ERROR << "Setting cross section to 1.0." <<endl;
          rh.setCrossSection(1.0);
        }
        #else
        log << Log::ERROR << "Warning: Cross section needed but no Generator created." <<endl;
        log << Log::ERROR << "Setting cross section to 1.0." <<endl;
        rh.setCrossSection(1.0);
        #endif
      }
      rh.finalize();
    }

    // Finalize generator
    //if (gen){
      ////Loader::destroyGen(gen);
      //Loader::finalize();
    //}

    // Write out the histogram tree into the registered file.
    rh.commitData();
  }

}
