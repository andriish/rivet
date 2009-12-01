// -*- C++ -*-
#include "Rivet/Run.hh"
#include "Rivet/AnalysisHandler.hh"
#include "HepMC/IO_GenEvent.h"
#include "Rivet/Projections/Beam.hh"
#include <limits>

namespace Rivet {


  Run::Run(AnalysisHandler& ah) 
    : _ah(ah), _xs(-1.0) 
  { }


  Run::~Run() { }


  Run& Run::setCrossSection(const double xs) {
    _xs = xs;
    return *this;
  }


  Run& Run::setListAnalyses(const bool dolist) {
    _listAnalyses = dolist;
    return *this;
  }


  bool Run::prepareFile(const std::string& evtfile) {
    if (evtfile == "-") {
      _io.reset(new HepMC::IO_GenEvent(std::cin));
    } else {
      // Ignore the HepMC::IO_GenEvent(filename, ios) constructor, since it's only available from HepMC 2.4
      _istr.reset(new std::fstream(evtfile.c_str(), std::ios::in));
      _io.reset(new HepMC::IO_GenEvent(*_istr));
    }
    if (_io->rdstate() != 0) {
      Log::getLog("Rivet.Run") << Log::ERROR << "Read error on file " << evtfile << endl;
      return false;
    }

    return true;
  }


  bool Run::processEvent(bool firstEvent) {
    // Fill event and check for a bad read state
    shared_ptr<GenEvent> evt;
    evt.reset(new GenEvent());
    if (_io->rdstate() != 0 || !_io->fill_next_event(evt.get()) ) {
      Log::getLog("Rivet.Run") << Log::DEBUG << "Read failed. End of file?" << endl;
      return false;
    }
 
    // Get beam details from first event, and ensure they match for all following events
    if (evt->particles_size() != 0) {
      const BeamPair beams = beamIds(*evt);
      const double sqrts = Rivet::sqrtS(*evt);
      Log::getLog("Rivet.Run") << Log::DEBUG << "Beams: "
                               << beams << " @ " << sqrts/GeV << " GeV" << endl;
      if (firstEvent) {
        _beams = beams;
        _sqrts = sqrts;
        Log::getLog("Rivet.Run") << Log::INFO << "First event beams: "
                                 << this->beams() << " @ " << this->sqrtS()/GeV << " GeV" << endl;
      } else {
        if (beams != _beams || !fuzzyEquals(sqrts, sqrtS())) {
          Log::getLog("Rivet.Run") << Log::ERROR << "Event beams mismatch: "
                                   << beams << " @ " << sqrts/GeV << " GeV" << " vs. first beams "
                                   << this->beams() << " @ " << this->sqrtS()/GeV << " GeV" << endl;
          return false;
        }
      }
    }

    // Set up system based on properties of first event
    if (firstEvent) {
      // If empty
      if (evt->particles_size() == 0) {
        Log::getLog("Rivet.Run") << Log::ERROR << "Empty first event." << endl;
        return false;
      }

      const size_t num_anas_requested = _ah.analysisNames().size();
      _ah.removeIncompatibleAnalyses(beams());
      if (num_anas_requested > 0 && _ah.analysisNames().size() == 0) {
        Log::getLog("Rivet.Run") << Log::ERROR
            << "All analyses were incompatible with the first event's beams\n"
            << "Exiting, since this probably isn't intentional!" << endl;
        return false;
      }
   
      if (_listAnalyses) {
        foreach (const std::string& ana, _ah.analysisNames()) {
          cout << ana << endl;
        }
      }

    }


    // Set cross-section if specified from command line
    if (_xs > 0.0) {
      if (firstEvent) {
        Log::getLog("Rivet.Run") << Log::DEBUG
                                 << "Setting user cross-section = " << _xs << " pb" << endl;
        _ah.setCrossSection(_xs);
      }
    }
    // Set cross-section if found in event
    #ifdef HEPMC_HAS_CROSS_SECTION
    else if (evt->cross_section()) {
      const double xs = evt->cross_section()->cross_section(); //< in pb
      Log::getLog("Rivet.Run") << Log::DEBUG
                               << "Setting cross-section = " << xs << " pb" << endl;
      _ah.setCrossSection(xs);
    }
    #endif
    // Complain about absence of cross-section if required!
    else {
      if (_ah.needCrossSection()) {
        Log::getLog("Rivet.Run") << Log::ERROR
            << "Total cross-section needed for at least one of the analyses. "
            << "Please set it (on the command line with '-x' if using the 'rivet' program)" << endl;
        return false;
      }
    }

    /// @todo If NOT first event, check that beams aren't changed
 
    // Analyze event
    _ah.analyze(*evt);
 
    return true;
  }


  bool Run::finalizeFile() {
    return true;
  }


  const BeamPair& Run::beams() const {
    return _beams;
  }


  double Run::sqrtS() const {
    return _sqrts;
  }


}
