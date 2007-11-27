// -*- C++ -*-
#include "AGILe/AGILe.hh"
#include "AGILe/Particle.hh"
#include "AGILe/Generator.hh"
using namespace AGILe;

#include "RivetGun/Configuration.hh"
#include "Rivet/Rivet.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Tools/Logging.hh"
using namespace Rivet;

#include "HepMC/IO_Ascii.h" // To become #include "HepMC/IO_GenEvent.h"
using namespace HepMC;


namespace Rivet {

  void generate(Generator* gen, Configuration& cfg, Log& log) {

    log << Log::DEBUG << "In generate method " << endl;

    // Seed random number generator
    gen->setSeed(cfg.rngSeed);

    bool needsCrossSection = false;
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
    } catch (runtime_error& e) {
      log << Log::ERROR << "Beam particle error: " << cfg.beam1 << " or " 
          << cfg.beam2 << " is invalid" << endl;
      throw runtime_error("Invalid beam particle");
    }

    log << Log::DEBUG << "Finished setting initial state " << endl;

    

    // Initialise Rivet
    AnalysisHandler rh(cfg.histoName, cfg.histoFormat); 
    if (cfg.runRivet) {
      for (set<string>::const_iterator ai = cfg.analyses.begin(); ai != cfg.analyses.end(); ++ai) {
        Analysis* a = AnalysisLoader::getAnalysis(*ai);
        log << Log::DEBUG << a->getName() << ": " << a->getBeams() << " "
            << (a->isCompatible(cfg.beam1, cfg.beam2) ? "MATCH" : "NO-MATCH") << endl;
        if (a->isCompatible(cfg.beam1, cfg.beam2)) {
          rh.addAnalysis(*ai);
          log << Log::INFO << "Using analysis " << a->getName() << endl;
        } else {
          log << Log::WARN << "Removing inappropriate analysis " << a->getName() << endl;
        }
        
        if (a->needsCrossSection()) { needsCrossSection = true; }
        
      }
      log << Log::DEBUG << "Initialising the Rivet analysis handler" << endl;
      rh.init();
      log << Log::DEBUG << "Rivet analysis handler initialised" << endl;
    }

    // Make a HepMC output strategy
    IO_Ascii* hepmcOut(0);
    if (cfg.writeHepMC) {
      /// @todo Use IO_GenEvent
      hepmcOut = new IO_Ascii(cfg.hepmcOutFile.c_str(), std::ios::out);
    }

    // Log the event number to a special logger
    Log& nevtlog = Log::getLog("RivetGun.NEvt");

    // Event loop
    log << Log::INFO << "Generating " << cfg.numEvents << " events."  << endl;
    HepMC::GenEvent myevent;
    for (unsigned int i = 0; i < cfg.numEvents; ++i) {    
      // Make the event
      gen->makeEvent(myevent);

      // Notify about event number
      Log::Level lev = Log::DEBUG;
      if (round((i+1)/100.0) == (i+1)/100.0) lev = Log::INFO;
      if (round((i+1)/1000.0) == (i+1)/1000.0) lev = Log::WARN;
      nevtlog << lev << "Event number " << i+1 << endl;

      // Run Rivet analyse s
      if (cfg.runRivet) rh.analyze(myevent);

      // Write out event to file
      if (cfg.writeHepMC) {
        hepmcOut->write_event(&myevent);
        myevent.print(cout);
      }
      /// @todo Clean-up
    }
    log << Log::INFO << "Finished!"  << endl;

    // Finalise Rivet and the generator
    gen->finalize();
    if (cfg.runRivet){
      if(needsCrossSection){
        rh.setCrossSection(gen->getCrossSection());
      }
      rh.finalize();
    }
  }

}
