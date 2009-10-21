// -*- C++ -*-
#include "Rivet/Run.hh"
#include "Rivet/AnalysisHandler.hh"
#include "HepMC/IO_GenEvent.h"
#include "Rivet/Projections/Beam.hh"
#include <limits>

namespace Rivet {


  Run::Run(AnalysisHandler& ah) : _ah(ah), _xs(-1.0),
    m_io(NULL), m_istr(NULL) {
  }
  
  
  Run::~Run() {
  }
  
  
  Run& Run::setCrossSection(const double& xs) {
    _xs = xs;
    return *this;
  }
  
  
  Run& Run::setListAnalyses(const bool& l) {
    _listAnalyses = l;
    return *this;
  }
  
  
  bool Run::prepareFile(const std::string& evtfile) {
    if (evtfile == "-") {
      m_io = new HepMC::IO_GenEvent(std::cin);
    } else {
      // Ignore the HepMC::IO_GenEvent(filename, ios) constructor, since only available from HepMC 2.4
      m_istr = new std::fstream(evtfile.c_str(), std::ios::in);
      m_io = new HepMC::IO_GenEvent(*m_istr);
    }
    if (m_io->rdstate() != 0) {
      Log::getLog("Rivet.Run") << Log::ERROR << "Read error on file " << evtfile << endl;
      return false;
    }
  
    return true;
  }
  
  bool Run::processEvent(bool firstEvent) {
    GenEvent* evt = new GenEvent();
    if (!m_io->fill_next_event(evt)) {
      delete evt;
      return false;
    }

    // Check for a bad read
    if (m_io->rdstate() != 0) {
      Log::getLog("Rivet.Run") << Log::DEBUG << "End of file?" << endl;
      delete evt;
      return false;
    }
    
    // Set up system based on properties of first event
    if (firstEvent) {
      // If empty
      if (evt->particles_size() == 0) {
        Log::getLog("Rivet.Run") << Log::ERROR << "Empty first event." << endl;
        delete evt;
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
        delete evt;
        return false;
      }
    }
    /// @todo If NOT first event, check that beams aren't changed
    
    // Analyze event and delete HepMC event object      
    _ah.analyze(*evt);
    delete evt;
    
    return true;
  }

  bool Run::finalizeFile() {
    // Final HepMC object clean-up
    delete m_io;
    if (m_istr) delete m_istr;
    
    return true;
  }

  
}
