// -*- C++ -*-
#include "AGILe/AGILe.hh"
#include "AGILe/Particle.hh"
#include "AGILe/Generator.hh"
using namespace AGILe;

#include "Rivet/Rivet.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/Configuration.hh"
using namespace Rivet;

#include "HepMC/IO_GenEvent.h"
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
        throw runtime_error("Couldn't read HepMC event file: " + cfg.hepmcInFile);
      }
    }
    
    // Log the event number to a special logger
    Log& nevtlog = Log::getLog("RivetGun.NEvt");
    
    // Event loop
    string rgverb = "Generating";
    if (cfg.readHepMC) rgverb = "Reading";
    log << Log::INFO << rgverb << " " << cfg.numEvents << " events." << endl;
    HepMC::GenEvent myevent;
    for (size_t i = 0; i < cfg.numEvents; ++i) {    
      // Make or load the event
      if (cfg.readHepMC) {
        if (hepmcIn->rdstate() != 0) {
          log << Log::ERROR << "Couldn't read next HepMC event from file: " << cfg.hepmcInFile << endl;
          break;
        }
        hepmcIn->fill_next_event(&myevent);
      } else {
        gen->makeEvent(myevent);
      }
      
      // Notify about event number
      Log::Level lev = Log::DEBUG;
      if (round((i+1)/100.0) == (i+1)/100.0) lev = Log::INFO;
      if (round((i+1)/1000.0) == (i+1)/1000.0) lev = Log::WARN;
      nevtlog << lev << "Event number " << i+1 << endl;
      
      // Run Rivet analyses
      if (cfg.runRivet) rh.analyze(myevent);
      
      // Write out event to file
      if (cfg.writeHepMC) {
        hepmcOut->write_event(&myevent);
        //myevent.print(cout);
      }
    }
    log << Log::INFO << "Finished!"  << endl;
    
    // Finalise Rivet and the generator
    gen->finalize();
    if (hepmcOut) delete hepmcOut;
    if (hepmcIn) delete hepmcIn;
    if (cfg.runRivet){
      if (needsCrossSection) {
        rh.setCrossSection(gen->getCrossSection());
      }
      rh.finalize();
    }
  }

}
