// -*- C++ -*-
#include "Rivet/Run.hh"
#include "Rivet/AnalysisHandler.hh"
#include "HepMC/IO_GenEvent.h"
#include "Rivet/Projections/Beam.hh"
#include <limits>

namespace Rivet {


  Run::Run(AnalysisHandler& ah) 
    : _ah(ah), _xs(-1.0), _sqrts(-1.0)
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


  // Fill event and check for a bad read state
  bool Run::readEvent() {
    /// @todo Clear rather than new the GenEvent object per-event?
    _evt.reset(new GenEvent());
    if (_io->rdstate() != 0 || !_io->fill_next_event(_evt.get()) ) {
      //Log::getLog("Rivet.Run") << Log::DEBUG << "Read failed. End of file?" << endl;
      return false;
    }
    return true;
  }


  bool Run::init(const std::string& evtfile) {
    // Set up HepMC input reader objects
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

    // Read first event to define run conditions
    bool ok = readEvent();
    if (!ok) return false;
    if (_evt->particles_size() == 0) {
      Log::getLog("Rivet.Run") << Log::ERROR << "Empty first event." << endl;
      return false;
    }

    // Set required beams for run based on first beams 
    const BeamPair beams = beamIds(*_evt);
    const double sqrts = Rivet::sqrtS(*_evt);
    _beams = beams;
    _sqrts = sqrts;
    Log::getLog("Rivet.Run") << Log::INFO << "First event beams: "
                             << this->beams() << " @ " << this->sqrtS()/GeV << " GeV" << endl;
    // Pass to analysis handler
    _ah.setBeams(_beams);
    _ah.setSqrtS(_sqrts);

    // Set cross-section from command line
    if (_xs >= 0.0) {
      Log::getLog("Rivet.Run") 
        << Log::DEBUG << "Setting user cross-section = " << _xs << " pb" << endl;
      _ah.setCrossSection(_xs);
    }

    // Check that analyses are beam-compatible
    const size_t num_anas_requested = _ah.analysisNames().size();
    _ah.removeIncompatibleAnalyses(beams);
    if (num_anas_requested > 0 && _ah.analysisNames().size() == 0) {
      Log::getLog("Rivet.Run") << Log::ERROR
                               << "All analyses were incompatible with the first event's beams\n"
                               << "Exiting, since this probably isn't intentional!" << endl;
      return false;
    }
   
    // List the chosen & compatible analyses if requested
    if (_listAnalyses) {
      foreach (const std::string& ana, _ah.analysisNames()) {
        cout << ana << endl;
      }
    }

    return true;
  }


  bool Run::processEvent() {
    // Ensure that beam details match those from first event
    const BeamPair beams = beamIds(*_evt);
    const double sqrts = Rivet::sqrtS(*_evt);
    if (beams != _beams || !fuzzyEquals(sqrts, sqrtS())) {
      Log::getLog("Rivet.Run") 
        << Log::ERROR << "Event beams mismatch: "
        << beams << " @ " << sqrts/GeV << " GeV" << " vs. first beams "
        << this->beams() << " @ " << this->sqrtS()/GeV << " GeV" << endl;
      return false;
    }

    // Set cross-section if found in event and not from command line
    #ifdef HEPMC_HAS_CROSS_SECTION
    if (_xs < 0.0 && _evt->cross_section()) {
      const double xs = _evt->cross_section()->cross_section(); //< in pb
      Log::getLog("Rivet.Run")
        << Log::DEBUG << "Setting cross-section = " << xs << " pb" << endl;
      _ah.setCrossSection(xs);
    }
    #endif
    // Complain about absence of cross-section if required!
    if (_ah.needCrossSection() && !_ah.hasCrossSection()) {
      Log::getLog("Rivet.Run") 
        << Log::ERROR
        << "Total cross-section needed for at least one of the analyses. "
        << "Please set it (on the command line with '-x' if using the 'rivet' program)" << endl;
      return false;
    }
     
    // Analyze event
    _ah.analyze(*_evt);
 
    return true;
  }


  bool Run::finalize() {
    _evt.reset();
    _istr.reset();
    _io.reset();
    return true;
  }


  const BeamPair& Run::beams() const {
    return _beams;
  }


  double Run::sqrtS() const {
    return _sqrts;
  }


}
