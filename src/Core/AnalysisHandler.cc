// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/ParticleName.hh"
#include "Rivet/Tools/BeamConstraint.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Beam.hh"
#include "YODA/WriterYODA.h"
#include <regex>

namespace {
    inline std::vector<std::string> split(const std::string& input, const std::string& regex) {
        // passing -1 as the submatch index parameter performs splitting
        std::regex re(regex);
        std::sregex_token_iterator
            first{input.begin(), input.end(), re, -1},
            last;
        return {first, last};
    }
}


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


  /// http://stackoverflow.com/questions/4654636/how-to-determine-if-a-string-is-a-number-with-c
  bool is_number(const std::string& s)
  {
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
  }
 
  /// Check if any of the weightnames is not a number 
  bool AnalysisHandler::haveNamedWeights() {
    bool dec=false;
    for (unsigned int i=0;i<_weightNames.size();++i) {
      string s = _weightNames[i];
      if (!is_number(s)) {
        dec=true;
        break;
      }
    }
    return dec;
  }


  void AnalysisHandler::init(const GenEvent& ge) {
    if (_initialised)
      throw UserError("AnalysisHandler::init has already been called: cannot re-initialize!");

    setRunBeams(Rivet::beams(ge));
    MSG_DEBUG("Initialising the analysis handler");
    _eventNumber = ge.event_number();

    setWeightNames(ge);
    if (haveNamedWeights())
        MSG_INFO("Using named weights");
    else
        MSG_INFO("NOT using named weights. Using first weight as nominal weight");

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
    /// reroute the print output to a stringstream and process
    /// The iteration is done over a map in hepmc2 so this is safe
    ostringstream stream;
    ge.weights().print(stream);  // Super lame, I know
    string str =  stream.str();

    std::regex re("(([^()]+))"); // Regex for stuff enclosed by parentheses ()
    for(std::sregex_iterator i = std::sregex_iterator(str.begin(), str.end(), re);
                          i != std::sregex_iterator(); ++i ) {
      std::smatch m = *i;
      vector<string> temp = split(m.str(), "[,]");
      if (temp.size() ==2) {
        MSG_DEBUG("Name of weight #" << _weightNames.size() << ": " << temp[0]);
        _weightNames.push_back(temp[0]);
      }
    }
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

        MSG_TRACE("Pushing analysis objects to persistent.");
        // if this is indeed a new event, push the temporary
        // histograms and reset
        for (const AnaHandle& a : _analyses) {
            for (const auto & ao : a->analysisObjects()) {
                ao.get().pushToPersistent(_subEventWeights);
            }
        }

        MSG_TRACE("Pushing eventCounter to persistent.");
        _eventCounter.pushToPersistent(_subEventWeights);

        _eventNumber = ge.event_number();

        MSG_DEBUG("Event #" << numEvents() << " weight sum = " << sumOfWeights());
        MSG_DEBUG("Event has " << _subEventWeights.size() << " sub events.");
        _subEventWeights.clear();
    }


    _eventCounter.newSubEvent();
    for (const AnaHandle& a : _analyses) {
        for (const auto & ao : a->analysisObjects()) {
            ao.get().newSubEvent();
        }
    }

    _subEventWeights.push_back(event.weights());
    MSG_DEBUG("Analyzing subevent #" << _subEventWeights.size() - 1 << ".");


    // Cross-section
    #ifdef HEPMC_HAS_CROSS_SECTION
    if (ge.cross_section()) {
      _xs = ge.cross_section()->cross_section();
      _xserr = ge.cross_section()->cross_section_error();
    }
    #endif

    MSG_TRACE("Filling event counter.");
    _eventCounter->fill();
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
          for (size_t iW = 0; iW < numWeights(); iW++) {
              for (MultiweightAOPtr& aoptr : a->_analysisobjects)
                  aoptr.setActiveWeightIdx(iW);

              try {
                  a->finalize();
              } catch (const Error& err) {
                  cerr << "Error in " << a->name() << "::finalize method: " << err.what() << endl;
                  exit(1);
              }
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


  vector<reference_wrapper<MultiweightAOPtr> > AnalysisHandler::getRivetAOs() const {
      vector<reference_wrapper<MultiweightAOPtr> > rtn;

      for (AnaHandle a : _analyses) {
          for (const auto & ao : a->analysisObjects()) {
              rtn.push_back(ao);
          }
      }

      // Should event counter be included here?
      rtn.push_back(_eventCounter);

      return rtn;
  }

  vector<YODA::AnalysisObjectPtr> AnalysisHandler::getYodaAOs() const {
      vector<YODA::AnalysisObjectPtr> rtn;

      const auto & raos = getRivetAOs();

      for (const auto & raoref : raos) {
          MultiweightAOPtr & rao = raoref.get();
          if (rao->path().find("/TMP/") != string::npos)
              continue;

          for (size_t iW = 0; iW < numWeights(); iW++) {
              rao.setActiveWeightIdx(iW);

              // add the weight name in brackets unless we recognize a
              // nominal weight
              if (_weightNames[iW] != "Weight" 
                  && _weightNames[iW] != "0"
                  && _weightNames[iW] != "Default")
                  rao->setPath(rao->path() + "[" + _weightNames[iW] + "]");

              rtn.push_back(rao.activeYODAPtr());
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
