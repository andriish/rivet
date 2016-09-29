// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/ParticleName.hh"
#include "Rivet/Tools/BeamConstraint.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Beam.hh"
#include "YODA/WriterYODA.h"


namespace Rivet {


  AnalysisHandler::AnalysisHandler(const string& runname)
    : _runname(runname),
      _eventCounter(0, Counter()), _xs(NAN),
      _initialised(false), _ignoreBeams(false)
  {}


  AnalysisHandler::~AnalysisHandler()
  {}


  Log& AnalysisHandler::getLog() const {
    return Log::getLog("Rivet.Analysis.Handler");
  }


  void AnalysisHandler::init(const GenEvent& ge) {
    if (_initialised)
      throw UserError("AnalysisHandler::init has already been called: cannot re-initialize!");

    setRunBeams(Rivet::beams(ge));
    MSG_DEBUG("Initialising the analysis handler");
    _eventNumber = ge.event_number();

    setWeightNames(ge);
    if (haveNamedWeights()) {
        MSG_INFO("Using named weights");
    }
    else {
        _weightNames.clear();
        _weightNames.push_back("0");
        MSG_INFO("NOT using named weights. Using first weight as nominal weight");
    }

    _numWeightTypes = _weightNames.size();
    _eventCounter = CounterPtr(_numWeightTypes, Counter("_EVTCOUNT"));

    // Check that analyses are beam-compatible, and remove those that aren't
    const size_t num_anas_requested = analysisNames().size();
    vector<string> anamestodelete;
    for (const AnaHandle a : _analyses) {
      if (!_ignoreBeams && !a->isCompatible(beams())) {
        //MSG_DEBUG(a->name() << " requires beams " << a->requiredBeams() << " @ " << a->requiredEnergies() << " GeV");
        anamestodelete.push_back(a->name());
      }
    }
    for (const string& aname : anamestodelete) {
      MSG_WARNING("Analysis '" << aname << "' is incompatible with the provided beams: removing");
      removeAnalysis(aname);
    }
    if (num_anas_requested > 0 && analysisNames().empty()) {
      cerr << "All analyses were incompatible with the first event's beams\n"
           << "Exiting, since this probably wasn't intentional!" << endl;
      exit(1);
    }

    // Warn if any analysis' status is not unblemished
    for (const AnaHandle a : analyses()) {
      if (toUpper(a->status()) == "PRELIMINARY") {
        MSG_WARNING("Analysis '" << a->name() << "' is preliminary: be careful, it may change and/or be renamed!");
      } else if (toUpper(a->status()) == "OBSOLETE") {
        MSG_WARNING("Analysis '" << a->name() << "' is obsolete: please update!");
      } else if (toUpper(a->status()).find("UNVALIDATED") != string::npos) {
        MSG_WARNING("Analysis '" << a->name() << "' is unvalidated: be careful, it may be broken!");
      }
    }

    // Initialize the remaining analyses
    for (AnaHandle a : _analyses) {
      MSG_DEBUG("Initialising analysis: " << a->name());
      try {
        // Allow projection registration in the init phase onwards
        a->_allowProjReg = true;
        a->init();
        //MSG_DEBUG("Checking consistency of analysis: " << a->name());
        //a->checkConsistency();
      } catch (const Error& err) {
        cerr << "Error in " << a->name() << "::init method: " << err.what() << endl;
        exit(1);
      }
      MSG_DEBUG("Done initialising analysis: " << a->name());
    }
    _initialised = true;
    MSG_DEBUG("Analysis handler initialised");
  }

  void AnalysisHandler::setWeightNames(const GenEvent& ge) {
      return;

      /// @todo
      /// this relies on boost and is removed for now
 
      /*
      /// reroute the print output to a stringstream and process
      /// The iteration is done over a map in hepmc2 so this is safe
      ostringstream stream;
      ge.weights().print(stream);  // Super lame, I know
      string str =  stream.str();

      vector<string> temp;
      split(temp, str, is_any_of("="));
      vector<string> tokens;
      split(tokens, temp[0], is_any_of(" "));

      /// Weights  come in tuples:  (StringName, doubleWeight) (AnotherName, anotherWeight) etc
      string wname;
      for (unsigned int i=0; i<tokens.size()-1;++i) {
        temp.clear();
        split(temp, tokens[i], is_any_of(","));
        if (temp.size()>0) {
          wname = temp[0];
          trim_left_if(wname, is_any_of("("));
          _weightNames.push_back(wname);
          /// Turn into debug message
          MSG_INFO("Name of weight #" << i << ": " << wname);
        }
      }
      */
  }


  void AnalysisHandler::analyze(const GenEvent& ge) {
    // Call init with event as template if not already initialised
    if (!_initialised) init(ge);
    assert(_initialised);

    // Ensure that beam details match those from the first event
    const PdgIdPair beams = Rivet::beamIds(ge);
    const double sqrts = Rivet::sqrtS(ge);
    if (!compatible(beams, _beams) || !fuzzyEquals(sqrts, sqrtS())) {
      cerr << "Event beams mismatch: "
           << PID::toBeamsString(beams) << " @ " << sqrts/GeV << " GeV" << " vs. first beams "
           << this->beams() << " @ " << this->sqrtS()/GeV << " GeV" << endl;
      exit(1);
    }

    // Create the Rivet event wrapper
    /// @todo Filter/normalize the event here
    Event event(ge);

    // won't happen for first event because _eventNumber is set in
    // init()
    if (_eventNumber != ge.event_number()) {
        /// @todo
        /// can we get away with not passing a matrix?

        // if this is indeed a new event, push the temporary
        // histograms and reset
        for (const AnaHandle& a : _analyses) {
            for (const shared_ptr<MultiweightAOPtr>& aoptr : a->analysisObjects()) {
                aoptr->pushToPersistent(_subEventWeights);
            }
        }

        _eventCounter.pushToPersistent(_subEventWeights);

        _eventNumber = ge.event_number();

        MSG_DEBUG("Event #" << numEvents() << " weight sum = " << sumOfWeights());
        MSG_DEBUG("Event has " << _subEventWeights.size() << " sub events.");
        _subEventWeights.clear();
    }


    for (const AnaHandle& a : _analyses) {
        for (const shared_ptr<MultiweightAOPtr>& aoptr : a->analysisObjects()) {
            aoptr->newSubEvent();
        }
    }

    _eventCounter.newSubEvent();

    _subEventWeights.push_back(event.weights());
    MSG_DEBUG("Analyzing subevent #" << _subEventWeights.size() - 1 << ".");


    // Cross-section
    #ifdef HEPMC_HAS_CROSS_SECTION
    if (ge.cross_section()) {
      _xs = ge.cross_section()->cross_section();
      _xserr = ge.cross_section()->cross_section_error();
    }
    #endif

    // Run the analyses
    for (AnaHandle a : _analyses) {
      MSG_TRACE("About to run analysis " << a->name());
      try {
        a->analyze(event);
      } catch (const Error& err) {
        cerr << "Error in " << a->name() << "::analyze method: " << err.what() << endl;
        exit(1);
      }
      MSG_TRACE("Finished running analysis " << a->name());
    }
  }


  void AnalysisHandler::analyze(const GenEvent* ge) {
    if (ge == NULL) {
      MSG_ERROR("AnalysisHandler received null pointer to GenEvent");
      //throw Error("AnalysisHandler received null pointer to GenEvent");
    }
    analyze(*ge);
  }


  void AnalysisHandler::finalize() {
    if (!_initialised) return;
    MSG_INFO("Finalising analyses");
    for (AnaHandle a : _analyses) {
      a->setCrossSection(_xs);
      try {
        a->finalize();
      } catch (const Error& err) {
        cerr << "Error in " << a->name() << "::finalize method: " << err.what() << endl;
        exit(1);
      }
    }

    // Print out number of events processed
    MSG_INFO("Processed " << numEvents() << " event" << (numEvents() == 1 ? "" : "s"));

    // // Delete analyses
    // MSG_DEBUG("Deleting analyses");
    // _analyses.clear();

    // Print out MCnet boilerplate
    cout << endl;
    cout << "The MCnet usage guidelines apply to Rivet: see http://www.montecarlonet.org/GUIDELINES" << endl;
    cout << "Please acknowledge plots made with Rivet analyses, and cite arXiv:1003.0694 (http://arxiv.org/abs/1003.0694)" << endl;
  }


  AnalysisHandler& AnalysisHandler::addAnalysis(const string& analysisname) {
    // Check for a duplicate analysis
    /// @todo Might we want to be able to run an analysis twice, with different params?
    ///       Requires avoiding histo tree clashes, i.e. storing the histos on the analysis objects.
    for (const AnaHandle& a : _analyses) {
      if (a->name() == analysisname) {
        MSG_WARNING("Analysis '" << analysisname << "' already registered: skipping duplicate");
        return *this;
      }
    }
    AnaHandle analysis( AnalysisLoader::getAnalysis(analysisname) );
    if (analysis.get() != 0) { // < Check for null analysis.
      MSG_DEBUG("Adding analysis '" << analysisname << "'");
      analysis->_analysishandler = this;
      _analyses.insert(analysis);
    } else {
      MSG_WARNING("Analysis '" << analysisname << "' not found.");
    }
    // MSG_WARNING(_analyses.size());
    // for (const AnaHandle& a : _analyses) MSG_WARNING(a->name());
    return *this;
  }


  AnalysisHandler& AnalysisHandler::removeAnalysis(const string& analysisname) {
    std::shared_ptr<Analysis> toremove;
    for (const AnaHandle a : _analyses) {
      if (a->name() == analysisname) {
        toremove = a;
        break;
      }
    }
    if (toremove.get() != 0) {
      MSG_DEBUG("Removing analysis '" << analysisname << "'");
      _analyses.erase(toremove);
    }
    return *this;
  }


  vector<shared_ptr<MultiweightAOPtr> > AnalysisHandler::getRivetAOs() const {
      vector<shared_ptr<MultiweightAOPtr> > rtn;

      for (AnaHandle a : _analyses) {
          for (const shared_ptr<MultiweightAOPtr>& aoptr : a->analysisObjects()) {
              rtn.push_back(aoptr);
          }
      }

      rtn.push_back(make_shared<CounterPtr>(_eventCounter));

      return rtn;
  }

  vector<YODA::AnalysisObjectPtr> AnalysisHandler::getYodaAOs() const {
      vector<YODA::AnalysisObjectPtr> rtn;

      vector<shared_ptr<MultiweightAOPtr> > raos =
          getRivetAOs();

      for (const shared_ptr<MultiweightAOPtr>& raoptr : raos) {
          if ((*raoptr)->path().find("/TMP/") != string::npos)
              continue;

          raoptr->setActiveWeightIdx(0);
          rtn.push_back(raoptr->activeYODAPtr());
          for (size_t iW = 1; iW < numWeights(); iW++) {
              raoptr->setActiveWeightIdx(iW);
              (*raoptr)->setPath((*raoptr)->path() + "[" + _weightNames[iW] + "]");
              rtn.push_back(raoptr->activeYODAPtr());
          }
      }

      YODA::Scatter1D::Points pts; pts.insert(YODA::Point1D(_xs, _xserr));
      rtn.push_back( make_shared<Scatter1D>(pts, "/_XSEC") );

      sort(rtn.begin(), rtn.end(),
           [](YODA::AnalysisObjectPtr a, YODA::AnalysisObjectPtr b) {
                return a->path() < b->path();
            }
          );

      return rtn;
  }

  vector<YODA::AnalysisObjectPtr> AnalysisHandler::getData() const {

      return getYodaAOs();
  }


  void AnalysisHandler::writeData(const string& filename) const {
    const vector<YODA::AnalysisObjectPtr> aos = getData();
    try {
      YODA::WriterYODA::write(filename, aos.begin(), aos.end());
    } catch ( YODA::WriteError ) {
      throw UserError("Unexpected error in writing file to: " + filename);
    }
  }


  string AnalysisHandler::runName() const { return _runname; }
  size_t AnalysisHandler::numEvents() const { return _eventCounter->numEntries(); }


  /*
   * why is this here?
  void AnalysisHandler::setSumOfWeights(const double& sum) {
    sumOfWeights() = sum;
  }
  */


  std::vector<std::string> AnalysisHandler::analysisNames() const {
    std::vector<std::string> rtn;
    for (AnaHandle a : _analyses) {
      rtn.push_back(a->name());
    }
    return rtn;
  }


  const AnaHandle AnalysisHandler::analysis(const std::string& analysisname) const {
    for (const AnaHandle a : analyses())
      if (a->name() == analysisname) return a;
    throw Error("No analysis named '" + analysisname + "' registered in AnalysisHandler");
  }


  AnalysisHandler& AnalysisHandler::addAnalyses(const std::vector<std::string>& analysisnames) {
    for (const string& aname : analysisnames) {
      //MSG_DEBUG("Adding analysis '" << aname << "'");
      addAnalysis(aname);
    }
    return *this;
  }


  AnalysisHandler& AnalysisHandler::removeAnalyses(const std::vector<std::string>& analysisnames) {
    for (const string& aname : analysisnames) {
      removeAnalysis(aname);
    }
    return *this;
  }


  bool AnalysisHandler::needCrossSection() const {
    bool rtn = false;
    for (const AnaHandle a : _analyses) {
      if (!rtn) rtn = a->needsCrossSection();
      if (rtn) break;
    }
    return rtn;
  }


  AnalysisHandler& AnalysisHandler::setCrossSection(double xs) {
    _xs = xs;
    return *this;
  }


  bool AnalysisHandler::hasCrossSection() const {
    return (!std::isnan(crossSection()));
  }


  AnalysisHandler& AnalysisHandler::addAnalysis(Analysis* analysis) {
    analysis->_analysishandler = this;
    _analyses.insert(AnaHandle(analysis));
    return *this;
  }


  PdgIdPair AnalysisHandler::beamIds() const {
    return Rivet::beamIds(beams());
  }


  double AnalysisHandler::sqrtS() const {
    return Rivet::sqrtS(beams());
  }

  void AnalysisHandler::setIgnoreBeams(bool ignore) {
    _ignoreBeams=ignore;
  }


}
