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
    std::istream* istr = 0;
    if (evtfile == "-") {
      io = new HepMC::IO_GenEvent(std::cin);
    } else {
      // Ignore the HepMC::IO_GenEvent(filename, ios) constructor, since only available from HepMC 2.4
      istr = new std::fstream(evtfile.c_str(), std::ios::in);
      io = new HepMC::IO_GenEvent(*istr);
    }
    if (io->rdstate() != 0) {
      Log::getLog("Rivet.Run") << Log::ERROR << "Read error on file " << evtfile << endl;
      return false;
    }
    
    GenEvent* evt = new GenEvent();
    while (io->fill_next_event(evt)) {

      // Check for a bad read
      if (io->rdstate() != 0) {
        Log::getLog("Rivet.Run") << Log::DEBUG << "End of file?" << endl;
        break;
      }

      // Set up system based on properties of first event
      if (_numEvents == 0) {
        // If empty
        if (evt->particles_size() == 0) {
          Log::getLog("Rivet.Run") << Log::ERROR << "Empty first event from " << evtfile << endl;
          return false;
        }
        size_t num_anas_requested = _ah.analysisNames().size();
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
              << "Please set it (on the command line with '-x' if using the 'rivet' program)" << endl;
          return false;
        }
      }
      /// @todo If NOT first event, check that beams aren't changed
      
      // Analyze event and delete HepMC event object      
      _ah.analyze(*evt);
      delete evt;

      // Increment, log, and check event number 
      ++_numEvents;
      logNEvt();
      if (_numEvents == _maxEvtNum) return false;
      
      /// @todo Can we just clear the event, to save on all this mallocing?
      evt = new GenEvent();
    }

    // Final HepMC object clean-up
    delete evt;
    delete io;
    delete istr;
    
    return true;
  }
  
  
  void Run::logNEvt() {
    std::stringstream ss;
    ss << "Event " << _numEvents;
    const string msg = ss.str();
    int lvl = Log::TRACE;
    if (_numEvents % 10000 == 0) lvl = Log::ERROR;
    else if (_numEvents % 1000 == 0) lvl = Log::WARN;
    else if (_numEvents % 100 == 0) lvl = Log::INFO;
    else if (_numEvents % 10 == 0) lvl = Log::DEBUG;
    Log::getLog("Rivet.Run") << lvl << msg << endl;
  }

  
}
