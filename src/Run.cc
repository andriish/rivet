// -*- C++ -*-
#include "Rivet/Run.hh"
#include "Rivet/AnalysisHandler.hh"
#include "HepMC/IO_GenEvent.h"
#include "Rivet/Projections/Beam.hh"
#include <limits>


namespace Rivet {

  Run::Run(AnalysisHandler& ah) : _ah(ah), _xs(-1.0),
    _maxEvtNum(std::numeric_limits<long>::max()), _numEvents(0) {
  }
  
  
  Run::~Run() {
    if (_maxEvtNum != std::numeric_limits<long>::max() && _numEvents < _maxEvtNum) {
      Log::getLog("Rivet.Run") << Log::WARN
          << "Sampled fewer events (" << _numEvents << ") than expected "
          << "(" << _maxEvtNum << ")" << endl;
    }
  }
  
  
  Run& Run::setCrossSection(const double& xs) {
    _xs = xs;
    return *this;
  }
  
  
  Run& Run::setMaxEvtNum(const int& n) {
    _maxEvtNum = n;
    return *this;
  }
  
  
  Run& Run::setListAnalyses(const bool& l) {
    _listAnalyses = l;
    return *this;
  }
  
  
  bool Run::processFile(const std::string& evtfile) {
    HepMC::IO_GenEvent* io = 0;
    if (evtfile == "-") {
      io = new HepMC::IO_GenEvent(std::cin);
    } else {
      io = new HepMC::IO_GenEvent(evtfile, std::ios::in);
    }
    if (io->rdstate() != 0) {
      Log::getLog("Rivet.Run") << Log::ERROR
          << "Read error on file " << evtfile << endl;
      return false;
    }
    
    GenEvent* evt = new GenEvent();
    while (io->fill_next_event(evt)) {
      if (_numEvents == 0) {
        int num_anas_requested = _ah.analysisNames().size();
        _ah.removeIncompatibleAnalyses(beamIds(*evt));
        if (num_anas_requested > 0 && _ah.analysisNames().size() == 0) {
          Log::getLog("Rivet.Run") << Log::ERROR
              << "All analyses were incompatible with the first event's beams"
              << "Exiting, since this probably isn't intentional!" << endl;
          delete evt;
          return false;
        }
        
        if (_listAnalyses) {
          foreach (const std::string& ana, _ah.analysisNames()) {
            cout<<ana<<endl;
          }
        }
      }
      
      if (_xs > 0.0) {
        _ah.setCrossSection(_xs);
      }
      #ifdef HEPMC_HAS_CROSS_SECTION
      else if (evt->cross_section()) {
        /// @todo Use xs error?
        const double xs = evt->cross_section()->cross_section(); //< in pb
        Log::getLog("Rivet.Run") << Log::DEBUG
            << "Setting cross-section = " << xs << " pb" << endl;
        _ah.setCrossSection(xs);
      }
      #endif
      else {
        if (_ah.needCrossSection()) {
          Log::getLog("Rivet.Run") << Log::ERROR
              << "Total cross-section needed for at least one of the analyses. "
              << "Please set it on the command line." << endl;
          return false;
        }
      }
      
      _ah.analyze(*evt);
      delete evt;
      
      ++_numEvents;
      logNEvt();
      
      if (_numEvents==_maxEvtNum) {
        return false;
      }
      
      evt = new GenEvent();
    }
    delete evt;
    delete io;
    
    return true;
  }
  
  
  void Run::logNEvt() {
    std::stringstream ss;
    ss << "Event " << _numEvents;
    if (_numEvents % 10 == 0)
      Log::getLog("Rivet.Run") << Log::DEBUG + 5 << ss.str() << endl;
    if (_numEvents % 100 == 0)
      Log::getLog("Rivet.Run") << Log::INFO << ss.str() << endl;
    if (_numEvents % 200 == 0)
      Log::getLog("Rivet.Run") << Log::INFO + 5 << ss.str() << endl;
    if (_numEvents % 500 == 0)
      Log::getLog("Rivet.Run") << Log::WARN << ss.str() << endl;
    if (_numEvents % 1000 == 0)
      Log::getLog("Rivet.Run") << Log::WARN + 5 << ss.str() << endl;
    if (_numEvents % 10000 == 0)
      Log::getLog("Rivet.Run") << Log::ERROR << ss.str() << endl;
    else
      Log::getLog("Rivet.Run") << Log::DEBUG << ss.str() << endl;
  }
  
}
