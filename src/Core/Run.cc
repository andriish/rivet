// -*- C++ -*-
#include "Rivet/Run.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Tools/RivetPaths.hh"
#include <limits>
#include <iostream>

using std::cout;
using std::endl;

namespace Rivet {

  Run::Run(AnalysisHandler& ah)
    : _ah(ah), _fileweight(1.0), _xs(NAN)
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
    if (!HepMCUtils::readEvent(_hepmcReader, _evt)){
      Log::getLog("Rivet.Run") << Log::DEBUG << "Read failed. End of file?" << endl;
      return false;
    }
    // Rescale event weights by file-level weight, if scaling is non-trivial
    if (_fileweight != 1.0) {
      for (size_t i = 0; i < (size_t) _evt->weights().size(); ++i) {
        _evt->weights()[i] *= _fileweight;
      }
    }
    return true;
  }



  bool Run::openFile(const std::string& evtfile, double weight) {
    // Set current weight-scaling member
    _fileweight = weight;

    // In case makeReader fails.
    std::string errormessage;

    // Use Rivet's own file format deduction (which uses the one in
    // HepMC3 if needed).
    _hepmcReader = HepMCUtils::makeReader(evtfile, _istr, &errormessage);

    // Check that it worked.
    if (_hepmcReader == nullptr) {
      Log::getLog("Rivet.Run")
        << Log::ERROR << "Read error in file '" << evtfile << "' "
        << errormessage << endl;
      return false;
    }
    return true;
  }


  bool Run::init(const std::string& evtfile, double weight) {
    if (!openFile(evtfile, weight)) return false;

    // Read first event to define run conditions
    bool ok = readEvent();
    if (!ok) return false;
    if(HepMCUtils::particles(_evt).size() == 0){
      Log::getLog("Rivet.Run") << Log::ERROR << "Empty first event." << endl;
      return false;
    }

    // Initialise AnalysisHandler with beam information from first event
    _ah.init(*_evt);

    // Set cross-section from command line
    if (!std::isnan(_xs)) {
      Log::getLog("Rivet.Run")
        << Log::DEBUG << "Setting user cross-section = " << _xs << " pb" << endl;

      _ah.setCrossSection(make_pair(_xs, 0.0));
    }

    // List the chosen & compatible analyses if requested
    if (_listAnalyses) {
      for (const std::string& ana : _ah.analysisNames()) {
        cout << ana << endl;
      }
    }

    return true;
  }


  bool Run::processEvent() {
    // Analyze event
    _ah.analyze(*_evt);

    return true;
  }


  bool Run::finalize() {
    _evt.reset();

    _ah.finalize();

    return true;
  }


}
