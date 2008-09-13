// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/HistoHandler.hh"
#include <algorithm>

namespace Rivet {


  // Initialize instance pointer to null.
  HistoHandler* HistoHandler::_instance = 0;


  HistoHandler* HistoHandler::create() {
    if (!_instance) {
      _instance = new HistoHandler();
      Log::getLog("Rivet.HistoHandler") << Log::TRACE << "Created new HistoHandler at " 
                                        << _instance << endl;
    }
    return _instance;
  }


  // Get a logger.
  Log& HistoHandler::getLog() const {
    return Log::getLog("Rivet.HistoHandler");
  }


  void HistoHandler::clear() {
    for (HistoHandles::iterator hh = _histos.begin(); hh != _histos.end(); ++hh) {
      getLog() << Log::TRACE << "Deleting projection at " << *hh << endl;
      delete *hh;
    }
    _histos.clear();
    _namedhistos.clear();
  }

  
  // Delete contained pointers.
  HistoHandler::~HistoHandler() {
    clear();
  }


  // Take a Projection, compare it to the others on record, and
  // return (by reference) an equivalent Projection which is guaranteed to be
  // the version that will be applied to an event.
  const Histo& HistoHandler::registerHisto(const Analysis& parent, 
                                           const AnalysisObject& ao, 
                                           const string& name) {
    getLog() << Log::TRACE << "Trying to register"
             << " analysis object " << &ao
             << " for parent " << &parent << "(" << parent.getName() << ")"
             << " with name '" << name << "'" << endl;
    

    /// @todo EVERYTHING!

  }



  const AnalysisObject& HistoHandler::getProjection(const Analysis& parent,
                                                    const string& name) const {
    getLog() << Log::TRACE << "Searching for child histo '" 
             << name << "' of " << &parent << endl;

    /// @todo EVERYTHING!

  }



}
