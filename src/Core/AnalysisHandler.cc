// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/ParticleName.hh"
#include "Rivet/Tools/BeamConstraint.hh"
#include "Rivet/Tools/RivetPaths.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Beam.hh"
#include "YODA/IO.h"
#include "YODA/WriterYODA.h"
#include <iostream>
#include <regex>
using namespace std;

namespace Rivet {


  AnalysisHandler::AnalysisHandler(const string& runname)
    : _runname(runname),
      _userxs{NAN, NAN},
      _initialised(false),
      _checkBeams(true),
      _skipMultiWeights(false),
      _matchWeightNames(""),
      _unmatchWeightNames(""),
      _nominalWeightName(""),
      _weightCap(0.),
      _NLOSmearing(0.), _defaultWeightIdx(0),
      _rivetDefaultWeightIdx(0), _dumpPeriod(0), _dumping(false)
  {  }


  AnalysisHandler::~AnalysisHandler() {
    static bool printed = false;
    // Print out MCnet boilerplate
    if (!printed && getLog().getLevel() <= 20) {
      cout << endl
           << "The MCnet usage guidelines apply to Rivet: see http://www.montecarlonet.org/GUIDELINES" << endl
           << "Please acknowledge Rivet in results made using it, and cite https://arxiv.org/abs/1912.05451" << endl;
      // << "https://arxiv.org/abs/1003.0694" << endl;
      printed = true;
    }
  }


  /// @todo Can we inline this?
  Log& AnalysisHandler::getLog() const {
    return Log::getLog("Rivet.AnalysisHandler");
  }


  /// http://stackoverflow.com/questions/4654636/how-to-determine-if-a-string-is-a-number-with-c
  namespace {
    bool is_number(const std::string& s) {
      std::string::const_iterator it = s.begin();
      while (it != s.end() && std::isdigit(*it)) ++it;
      return !s.empty() && it == s.end();
    }
  }

  /// Check if any of the weight names is not a number
  bool AnalysisHandler::haveNamedWeights() const {
    bool dec = false;
    for (size_t i = 0; i <_weightNames.size(); ++i) {
      const string& s = _weightNames[i];
      if (!is_number(s)) {
        dec = true;
        break;
      }
    }
    return dec;
  }


  void AnalysisHandler::init(const GenEvent& ge) {
    if (_initialised)
      throw UserError("AnalysisHandler::init has already been called: cannot re-initialize!");

    // Set the Run's beams based on this first event
    /// @todo Improve this const/ptr GenEvent mess!
    const Event e(const_cast<GenEvent&>(ge));
    setRunBeams(beams(e));
    MSG_DEBUG("Initialising the analysis handler");
    _eventNumber = ge.event_number();

    // Assemble the weight streams to be used
    setWeightNames(ge);
    if (_skipMultiWeights) {
      MSG_INFO("Only using nominal weight: variation weights will be ignored");
    } else if (haveNamedWeights()) {
      MSG_DEBUG("Using named weights");
    } else {
      MSG_WARNING("NOT using named weights: assuming first weight is nominal");
    }

    // Create the multi-weighted event counter
    _eventCounter = CounterPtr(weightNames(), Counter("_EVTCOUNT"));

    // Set the cross section based on what is reported by the init-event, else zero
    if (ge.cross_section()) {
      setCrossSection(HepMCUtils::crossSection(ge));
    } else {
      MSG_DEBUG("No cross-section detected in first event: setting default to 0 pb");
      setCrossSection({0.0, 0.0});
    }

    // Check that analyses are beam-compatible, and note those that aren't
    const size_t num_anas_requested = analysisNames().size();
    vector<string> anamestodelete;
    for (const AnaHandle& a : analyses()) {
      if (_checkBeams && !a->compatibleWithRun()) anamestodelete += a->name();
    }
    // Remove incompatible analyses from the run
    for (const string& aname : anamestodelete) {
      MSG_WARNING("Analysis '" << aname << "' is incompatible with the provided beams: removing");
      removeAnalysis(aname);
    }
    if (num_anas_requested > 0 && analysisNames().empty()) {
      MSG_ERROR("All analyses were incompatible with the first event's beams\n"
                << "Exiting, since this probably wasn't intentional!");
      exit(1);
    }

    // Warn if any analysis' status is not unblemished
    for (const AnaHandle& a : analyses()) {
      if ( a->info().preliminary() ) {
        MSG_WARNING("Analysis '" << a->name() << "' is preliminary: be careful, it may change and/or be renamed!");
      } else if ( a->info().obsolete() ) {
        MSG_WARNING("Analysis '" << a->name() << "' is obsolete: please update!");
      } else if (( a->info().unvalidated() ) ) {
        MSG_WARNING("Analysis '" << a->name() << "' is unvalidated: be careful, it may be broken!");
      }
    }

    // Initialize the remaining analyses
    _stage = Stage::INIT;
    for (AnaHandle a : analyses()) {
      MSG_DEBUG("Initialising analysis: " << a->name());
      try {
        // Allow projection registration in the init phase onwards
        a->_allowProjReg = true;
        a->setProjectionHandler(_projHandler);
        a->init();
        a->syncDeclQueue();
        //MSG_DEBUG("Checking consistency of analysis: " << a->name());
        //a->checkConsistency();
      } catch (const Error& err) {
        cerr << "Error in " << a->name() << "::init method: " << err.what() << endl;
        exit(1);
      }
      MSG_DEBUG("Done initialising analysis: " << a->name());
    }
    _stage = Stage::OTHER;
    _initialised = true;
    MSG_DEBUG("Analysis handler initialised");
  }


  // void AnalysisHandler::init(GenEvent* ge) {
  //   if (ge == nullptr) {
  //     MSG_ERROR("AnalysisHandler received null pointer to GenEvent");
  //     //throw Error("AnalysisHandler received null pointer to GenEvent");
  //   }
  //   init(*ge);
  // }


  void AnalysisHandler::setWeightNames(const GenEvent& ge) {
    _weightNames = HepMCUtils::weightNames(ge);

    // If there are no weights, add a nominal one
    if (_weightNames.empty()) {
      _weightNames.push_back("");
      _rivetDefaultWeightIdx = _defaultWeightIdx = 0;
      _weightIndices = { 0 };
      return;
    }

    // Find default weights, starting with the chosen or preferred name (default = "")
    size_t nDefaults = 0;
    _weightIndices.clear();
    for (size_t i = 0, N = _weightNames.size(); i < N; ++i) {
      _weightIndices.push_back(i);
      if (_weightNames[i] == _nominalWeightName) {
        if (nDefaults == 0) {
          _weightNames[i] = "";
          _rivetDefaultWeightIdx = _defaultWeightIdx = i;
        }
        nDefaults += 1;
      }
    }

    // If there are no weights with the preferred name, look for acceptable alternatives
    if (nDefaults == 0) {
      for (size_t i = 0, N = _weightNames.size(); i < N; ++i) {
        const string W = toUpper(_weightNames[i]);
        if (W == "WEIGHT" || W == "0" || W == "DEFAULT" || W == "NOMINAL") {
          if (nDefaults == 0) {
            _weightNames[i] = "";
            _rivetDefaultWeightIdx = _defaultWeightIdx = i;
          }
          nDefaults += 1;
        }
      }
    }

    // Warn user that no nominal weight could be identified
    if (nDefaults == 0) {
      MSG_WARNING("Could not identify nominal weight. Will continue assuming variations-only run.");
      // Note quoting for clarity, given the indents:
      MSG_WARNING("Candidate weight names:\n    '" << join(_weightNames, "'\n    '") << "'");
    }
    // Warn if multiple weight names were acceptable alternatives
    if (nDefaults > 1) {
      MSG_WARNING("Found more than " << nDefaults << " default weight candidates. Will use: " << _weightNames[_defaultWeightIdx]);
    }

    // Apply behaviours for only using the nominal weight, or all weights
    if (_skipMultiWeights)  {

      // If running in single-weight mode, remove all bar the nominal weight
      _weightIndices = { _defaultWeightIdx };
      _weightNames = { _weightNames[_defaultWeightIdx] };
      _rivetDefaultWeightIdx = 0;

    } else {

      // Check if weight name matches a supplied string/regex and filter to select those only
      if (_matchWeightNames != "") {
        MSG_DEBUG("Select weight names that match pattern \"" << _matchWeightNames << "\"");
        // Compile regex from each string in the comma-separated list
        vector<std::regex> patterns;
        for (const string& pattern : split(_matchWeightNames, ",")) {
          patterns.push_back( std::regex(pattern) );
        }
        // Check which weights match supplied weight-name pattern
        vector<string> selected_subset; vector<size_t> selected_indices;
        for (size_t i = 0, N = _weightNames.size(); i < N; ++i) {
          if (_weightIndices[i] == _defaultWeightIdx) {
            // The default weight cannot be "unselected"
            _rivetDefaultWeightIdx = selected_indices.size();
            selected_indices.push_back(_weightIndices[i]);
            selected_subset.push_back(_weightNames[i]);
            MSG_DEBUG("Selected nominal weight: \"" << _weightNames[i] << "\"");
            continue;
          }
          for (const std::regex& re : patterns) {
            if ( std::regex_match(_weightNames[i], re) ) {
              selected_indices.push_back(_weightIndices[i]);
              selected_subset.push_back(_weightNames[i]);
              MSG_DEBUG("Selected variation weight: \"" << _weightNames[i] << "\"");
              break;
            }
          }
        }
        _weightNames = selected_subset;
        _weightIndices = selected_indices;
      }

      // Check if the remaining weight names match supplied string/regexes and *de*select accordingly
      vector<std::regex> patterns = { std::regex("AUX"), std::regex("DEBUG") };
      if (_unmatchWeightNames != "") {
        MSG_DEBUG("Deselect weight names that match pattern \"" << _unmatchWeightNames << "\"");
        // Compile regex from each string in the comma-separated list
        for (const string& pattern : split(_unmatchWeightNames, ",")) {
          patterns.push_back( std::regex(pattern) );
        }
      }
      // Check which weights match supplied weight-name pattern
      vector<string> selected_subset; vector<size_t> selected_indices;
      for (size_t i = 0, N = _weightNames.size(); i < N; ++i) {
        if (_weightIndices[i] == _defaultWeightIdx) {
          // The default weight cannot be vetoed
          _rivetDefaultWeightIdx = selected_indices.size();
          selected_indices.push_back(_weightIndices[i]);
          selected_subset.push_back(_weightNames[i]);
          MSG_DEBUG("Selected nominal weight: " << _weightNames[i]);
          continue;
        }
        bool skip = false;
        for (const std::regex& re : patterns) {
          if ( std::regex_match(_weightNames[i], re) ) { skip = true; break; }
        }
        if (skip) continue;
        selected_indices.push_back(_weightIndices[i]);
        selected_subset.push_back(_weightNames[i]);
        MSG_DEBUG("Selected variation weight: " << _weightNames[i]);
      }
      _weightNames = selected_subset;
      _weightIndices = selected_indices;

    }

    // Done (de-)selecting weights: show useful debug messages
    MSG_DEBUG("Default weight name: \"" <<  _weightNames[_rivetDefaultWeightIdx] << "\"");
    MSG_DEBUG("Default weight index (Rivet): " << _rivetDefaultWeightIdx);
    MSG_DEBUG("Default weight index (overall): " << _defaultWeightIdx);
  }


  // bool AnalysisHandler::consistentWithRun(Event& event) {
  //   const PdgIdPair beamids = beamIDs(event);
  //   const double sqrts = sqrtS(event);
  //   return compatibleBeams(beamids, runBeamIDs()) && compatibleBeamEnergy(sqrts, runSqrtS());
  // }


  void AnalysisHandler::analyze(GenEvent& ge) {
    // Call init with event as template if not already initialised
    if (!_initialised) init(ge);
    assert(_initialised);

    // Create the Rivet event wrapper
    //bool strip = ( getEnvParam("RIVET_STRIP_HEPMC", string("NOOOO") ) != "NOOOO" );
    const Event event(ge, _weightIndices);

    // Ensure that beam details match those from the first event (if we're checking beams)
    if (_checkBeams) {
      const ParticlePair evtbeams = beams(event);
      MSG_DEBUG("Event beams = " << evtbeams);
      if (evtbeams.first.pid() == PID::ANY && evtbeams.second.pid() == PID::ANY) {
        MSG_ERROR("No event beams found: please fix the events, or run with beam-checking disabled");
        exit(1);
      }
      if (!compatibleBeams(evtbeams, runBeams())) {
        MSG_ERROR("Event beams mismatch with run: "
                  << PID::toBeamsString(beamIDs(event)) << " @ " << sqrtS(event)/GeV << " GeV" << " vs. expected "
                  << this->runBeams() << " @ " << this->runSqrtS()/GeV << " GeV");
        exit(1);
      }
    }

    // Set the cross section based on what is reported by this event
    if (ge.cross_section()) setCrossSection(HepMCUtils::crossSection(ge));

    // If the event number has changed, sync the sub-event analysis objects to persistent
    // NB. Won't happen for first event because _eventNumber is set in init()
    /// @todo Need to be able to turn this off, in the case of slightly malformed events without event numbers
    if (_eventNumber != ge.event_number()) {
      pushToPersistent();
      _eventNumber = ge.event_number();
    }

    // Make a new sub-event: affects every analysis object
    MSG_TRACE("Starting new sub-event");
    _eventCounter.get()->newSubEvent();
    for (const AnaHandle& a : analyses()) {
      for (auto ao : a->analysisObjects()) {
        ao.get()->newSubEvent();
      }
    }

    // Optionally cap large weights, to avoid spikes
    _subEventWeights.push_back(event.weights());
    if (_weightCap != 0.) {
      MSG_DEBUG("Implementing weight cap using a maximum |weight| = " << _weightCap << " for latest subevent");
      size_t lastSub = _subEventWeights.size() - 1;
      for (size_t i = 0; i < _subEventWeights[lastSub].size(); ++i) {
        if (abs(_subEventWeights[lastSub][i]) > _weightCap) {
          _subEventWeights[lastSub][i] = sign(_subEventWeights[lastSub][i]) * _weightCap;
        }
      }
    }

    // Run the analyses
    if (_subEventWeights.size() > 1) {
      MSG_TRACE("Analyzing subevent #" << _subEventWeights.size());
    }
    _eventCounter->fill();
    for (AnaHandle a : analyses()) {
      MSG_TRACE("About to run analysis " << a->name());
      try {
        a->analyze(event);
      } catch (const Error& err) {
        cerr << "Error in " << a->name() << "::analyze method: " << err.what() << endl;
        exit(1);
      }
      MSG_TRACE("Finished running analysis " << a->name());
    }

    // Dump current final histograms
    if ( _dumpPeriod > 0 && numEvents() > 0 && numEvents() % _dumpPeriod == 0 ) {
      MSG_DEBUG("Dumping intermediate results to " << _dumpFile << ".");
      _dumping = numEvents()/_dumpPeriod;
      finalize();
      writeData(_dumpFile);
      _dumping = 0;
    }

  }


  void AnalysisHandler::analyze(GenEvent* ge) {
    if (ge == nullptr) {
      MSG_ERROR("AnalysisHandler received null pointer to GenEvent");
      //throw Error("AnalysisHandler received null pointer to GenEvent");
    }
    analyze(*ge);
  }


  void AnalysisHandler::pushToPersistent() {
    if ( _subEventWeights.empty() ) return;
    MSG_TRACE("Count at " << __FILE__ <<": "<<__LINE__ << " is " << getEventCounter());
    MSG_TRACE("AnalysisHandler::analyze(): Pushing _eventCounter to persistent.");
    _eventCounter.get()->pushToPersistent(_subEventWeights);
    for (const AnaHandle& a : analyses()) {
      for (auto ao : a->analysisObjects()) {
        MSG_TRACE("AnalysisHandler::analyze(): Pushing " << a->name()
                  << "'s " << ao->name() << " to persistent.");
        ao.get()->pushToPersistent(_subEventWeights, _NLOSmearing);
      }
      MSG_TRACE("AnalysisHandler::analyze(): finished pushing "
                << a->name() << "'s objects to persistent.");
    }
    _subEventWeights.clear();
  }


  void AnalysisHandler::finalize() {
    if (!_initialised) return;
    MSG_DEBUG("Finalising analyses");

    _stage = Stage::FINALIZE;

    MSG_TRACE("Count at " << __FILE__ <<": "<<__LINE__ << " is " << getEventCounter());

    // First push all analyses' objects to persistent and final
    MSG_TRACE("AnalysisHandler::finalize(): Pushing analysis objects to persistent.");
    pushToPersistent();

    // Warn if no cross-section was set
    if (!nominalCrossSection()) {
      MSG_WARNING("Null nominal cross-section: setting to 10^-10 pb to allow rescaling");
      setCrossSection(1.0e-10, 0.0);
      // _xs.get()->setActiveWeightIdx(_rivetDefaultWeightIdx);
      // _xs->point(0).setX(1.0, 1.0);
      // _xs.get()->unsetActiveWeight();
    }

    // Copy all histos to finalize versions.
    _eventCounter.get()->pushToFinal();
    _xs.get()->pushToFinal();
    for (const AnaHandle& a : analyses()) {
      for (auto ao : a->analysisObjects()) {
        ao.get()->pushToFinal();
      }
    }

    for (AnaHandle a : analyses()) {
      if ( _dumping && !a->info().reentrant() )  {
        if ( _dumping == 1 )
          MSG_DEBUG("Skipping finalize in periodic dump of " << a->name() << " as it is not declared re-entrant.");
        continue;
      }
      for (size_t iW = 0; iW < numWeights(); iW++) {
        _eventCounter.get()->setActiveFinalWeightIdx(iW);

        MSG_TRACE("Count at " << __FILE__ <<": "<<__LINE__ << " is " << getEventCounter());
        _xs.get()->setActiveFinalWeightIdx(iW);
        for (auto ao : a->analysisObjects()) {
          ao.get()->setActiveFinalWeightIdx(iW);
        }
        try {
          MSG_TRACE("running " << a->name() << "::finalize() for weight " << iW << ".");
          a->finalize();
        } catch (const Error& err) {
          cerr << "Error in " << a->name() << "::finalize method: " << err.what() << '\n';
          exit(1);
        }
      }
    }

    // Print out number of events processed
    MSG_TRACE("Count at " << __FILE__ <<": "<<__LINE__ << " is " << getEventCounter());
    _eventCounter.get()->setActiveFinalWeightIdx(defaultWeightIndex());
    MSG_TRACE("Count at " << __FILE__ <<": "<<__LINE__ << " is " << getEventCounter());
    _xs.get()->setActiveFinalWeightIdx(defaultWeightIndex());
    MSG_TRACE("Count at " << __FILE__ <<": "<<__LINE__ << " is " << getEventCounter());
    if (!_dumping) {
      const int nevts = numEvents();
      MSG_DEBUG("Processed " << nevts << " event" << (nevts != 1 ? "s" : ""));
    }

    _stage = Stage::OTHER;

    MSG_TRACE("Count at " << __FILE__ <<": "<<__LINE__ << " is " << getEventCounter());

  }


  AnalysisHandler& AnalysisHandler::addAnalysis(const string& analysisname, std::map<string, string> pars) {
     // Make an option handle.
    std::string parHandle = "";
    for (map<string, string>::iterator par = pars.begin(); par != pars.end(); ++par) {
      parHandle +=":";
      parHandle += par->first + "=" + par->second;
    }
    return addAnalysis(analysisname + parHandle);
  }


  AnalysisHandler& AnalysisHandler::addAnalysis(const string& analysisname) {
    // Check for a duplicate analysis
    /// @todo Might we want to be able to run an analysis twice, with different params?
    ///       Requires avoiding histo tree clashes, i.e. storing the histos on the analysis objects.
    string ananame = analysisname;
    vector<string> anaopt = split(analysisname, ":");
    if ( anaopt.size() > 1 ) ananame = anaopt[0];
    AnaHandle analysis( AnalysisLoader::getAnalysis(ananame) );
    if (analysis.get() != 0) { // < Check for null analysis.
      MSG_DEBUG("Adding analysis '" << analysisname << "'");
      map<string,string> opts;
      for ( int i = 1, N = anaopt.size(); i < N; ++i ) {
        vector<string> opt = split(anaopt[i], "=");
        if ( opt.size() != 2 ) {
          MSG_WARNING("Error in option specification. Skipping analysis " << analysisname);
          return *this;
        }
        if ( !analysis->info().validOption(opt[0], opt[1]) )
          MSG_WARNING("Setting the option '" << opt[0] << "' to '"
                      << opt[1] << "' for " << analysisname
                      << " has not been declared in the info file "
                      << " and may be ignored in the analysis.");
        opts[opt[0]] = opt[1];
      }
      for ( auto opt: opts) {
        analysis->_options[opt.first] = opt.second;
        analysis->_optstring += ":" + opt.first + "=" + opt.second;
      }
      for (const AnaHandle& a : analyses()) {
        if (a->name() == analysis->name() ) {
          MSG_WARNING("Analysis '" << analysisname << "' already registered: skipping duplicate");
          return *this;
        }
      }
      analysis->_analysishandler = this;
      _analyses[analysisname] = analysis;
    } else {
      MSG_WARNING("Analysis '" << analysisname << "' not found.");
    }
    // MSG_WARNING(_analyses.size());
    // for (const AnaHandle& a : _analyses) MSG_WARNING(a->name());
    return *this;
  }


  AnalysisHandler& AnalysisHandler::removeAnalysis(const string& analysisname) {
    MSG_DEBUG("Removing analysis '" << analysisname << "'");
    if (_analyses.find(analysisname) != _analyses.end()) _analyses.erase(analysisname);
    // }
    return *this;
  }


  void AnalysisHandler::stripOptions(YODA::AnalysisObjectPtr ao,
                                     const vector<string> & delopts) const {
    string path = ao->path();
    string ananame = split(path, "/")[0];
    vector<string> anaopts = split(ananame, ":");
    for ( int i = 1, N = anaopts.size(); i < N; ++i )
      for ( auto opt : delopts )
        if ( opt == "*" || anaopts[i].find(opt + "=") == 0 )
          path.replace(path.find(":" + anaopts[i]), (":" + anaopts[i]).length(), "");
    ao->setPath(path);
  }


  void AnalysisHandler::mergeYodasFromFiles(const vector<string> &aofiles,
                                            const vector<string> &delopts,
                                            const vector<string> &addopts,
                                            const vector<string> &matches,
                                            const vector<string> &unmatches,
                                            bool equiv) {

    // Parse option adding.
    vector<string> optAnas;
    vector<string> optKeys;
    vector<string> optVals;
    for (string addopt : addopts) {
      size_t pos1 = addopt.find(":");
      size_t pos2 = addopt.find("=");
      if (pos1 == string::npos || pos2 == string::npos || pos2 < pos1) {
        MSG_WARNING("Malformed analysis option: "+addopt+". Format as ANA:OPT=VAL");
        continue;
      }
      optAnas.push_back(addopt.substr(0, pos1));
      optKeys.push_back(addopt.substr(pos1 +1, pos2 - pos1 - 1));
      optVals.push_back(addopt.substr(pos2 +1 , addopt.size() - pos2 - 1));
    }

    // Go through all files and collect information
    /// @todo Move this to the script interface, with the API working in terms
    ///   of <real_filename,weight> pairs rather than decoding a CLI convention in C++
    bool overwrite_xsec = false;
    size_t nfiles = 0, nfilestot = aofiles.size();
    map<string, YODA::AnalysisObjectPtr> allaos;
    map<string, pair<double,double> > allxsecs;
    for (string file : aofiles) {
      ++nfiles;
      std::cout << "Merging data file " << file << " [" << nfiles << "/" << nfilestot << "]\r";
      std::cout.flush();
      MSG_DEBUG("Reading in data from " << file);

      // Check for user-supplied scaling, assign 1 otherwise
      /// @todo
      size_t colonpos = file.rfind(":");
      double fileweight = 1.0;
      if (colonpos != string::npos) {
        string suffix = file.substr(colonpos+1);
        try {
          if (suffix.at(0) == '=') {
            // case I: file.yoda:=1.23
            //-> set cross-section to 1.23
            overwrite_xsec = true;
            suffix = suffix.substr(1);
          }
          else if (suffix.at(0) == 'x') {
            // case II: file.yoda:x1.23
            // (same as file.yoda:1.23)
            //-> multiply cross-section with 1.23
            suffix = suffix.substr(1);
          }
          fileweight = std::stod(suffix);
          file = file.substr(0, colonpos);
        } catch (...) {
          throw UserError("Unexpected error in processing argument " + file + " with file:scale format");
        }
      }

      // try to read the file and build path-AO map
      // @todo move this map construction into YODA?
      vector<YODA::AnalysisObject*> aos_raw;
      map<string,YODA::AnalysisObject*> raw_map;
      try {
        YODA::read(file, aos_raw);
        for (YODA::AnalysisObject* aor : aos_raw) {
          const string& aopath = aor->path();
          bool skip = false;
          if (aopath != "") {
            if (matches.size()) {
              skip = !std::any_of(matches.begin(), matches.end(), [&](const string &exp){
                                  return std::regex_match(aopath, std::regex(exp));} );
            }
            if (unmatches.size()) {
              skip |= std::any_of(unmatches.begin(), unmatches.end(), [&](const string &exp){
                                 return std::regex_match(aopath, std::regex(exp));} );
            }
          }
          if (skip)  continue;
          raw_map[aopath] = aor;
        }
      }
      catch (...) { //< YODA::ReadError&
        throw UserError("Unexpected error in reading file: " + file);
      }
      if (raw_map.empty()) {
        MSG_WARNING("No AOs selected from file: " << file);
        continue;
      }

      // merge AOs from current file into "allaos"
      mergeAOS(allaos, raw_map, allxsecs, delopts, optAnas, optKeys, optVals,
                                equiv, overwrite_xsec, fileweight);

    } // loop over all input files ends
    std::cout << std::endl;

    MSG_INFO("Rerunning finalize ...");

    // initialise analyses and load merged AOs back into memory
    setupReentrantRun(allaos, allxsecs, equiv);

    // Finally we just have to finalize all analyses, leaving to the
    // controlling program to write it out to some YODA file.
    finalize();
  }


  void AnalysisHandler::mergeAOS(map<string, YODA::AnalysisObjectPtr> &allaos,
                                 map<string, YODA::AnalysisObject*> &newaos, 
                                 map<string, pair<double, double>> &allxsecs,
                                 const vector<string> &delopts,
                                 const vector<string> &optAnas,
                                 const vector<string> &optKeys,
                                 const vector<string> &optVals,
                                 const bool equiv, 
                                 const bool overwrite_xsec,
                                 const double fileweight) {


    for (const auto & ao : newaos){
      std::cout << ao.first << ", " << ao.second << " (" << ao.second->path() << ")" << std::endl;
    }
                                   
    std::cout << __FILE__ << ": " << __LINE__ << std::endl;
    map<string, double> scales;
    for (auto& [aopath, aor] : newaos) { 
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      YODA::AnalysisObjectPtr ao(aor);
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      //AOPath path(ao->path());
      AOPath path(aopath);
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      if ( !path ) {
        throw UserError("Invalid path name in new AO set!");
      }
      // skip everything that isn't pre-finalize
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      if ( !path.isRaw() ){
        std::cout << __FILE__ << ": " << __LINE__ << std::endl;
        continue;
      } 

      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      MSG_DEBUG(" " << ao->path());

      const string& wname = path.weightComponent();
      if ( scales.find(wname) == scales.end() ) {
        std::cout << __FILE__ << ": " << __LINE__ << std::endl;
        scales[wname] = 1.0;
        // get the sum of weights and number of entries for the current weight
        double evts = 0, sumw = 1;
        auto ec_it = newaos.find("/RAW/_EVTCOUNT" + wname);
        if ( ec_it != newaos.end() ) {
          YODA::Counter* cPtr = static_cast<YODA::Counter*>(ec_it->second);
          evts = cPtr->numEntries();
          sumw = cPtr->sumW()? cPtr->sumW() : 1;
        }
        else if (!equiv) {
          throw UserError("Missing event counter, needed for non-equivalent merging!");
        }
        std::cout << __FILE__ << ": " << __LINE__ << std::endl;
        // in stacking mode: add up all the cross sections
        // in equivalent mode: weight the cross-sections
        // estimates by the corresponding number of entries
        const string xspath = "/RAW/_XSEC" + wname;
        auto xs_it = newaos.find(xspath);
        if ( xs_it != newaos.end() ) {
          YODA::Scatter1D* xsec = static_cast<YODA::Scatter1D*>(xs_it->second);
          if (overwrite_xsec) {
            MSG_DEBUG("Set user-supplied weight: " << fileweight);
            xsec->point(0).setX(fileweight);
          }
          else {
            MSG_DEBUG("Multiply user-supplied weight: " << fileweight);
            xsec->scaleX(fileweight);
          }
          // get iterator to the existing (or newly created) key-value pair
          auto xit = allxsecs.insert( make_pair(xspath, make_pair(0,0)) ).first;
          // update cross-sections, possibly weighted by number of entries
          xit->second.first  += (equiv? evts : 1.0) * xsec->point(0).x();
          xit->second.second += (equiv? sqr(evts) : 1.0) * sqr(xsec->point(0).xErrAvg());
          // only in stacking mode: multiply each AO by cross-section / sumW
          if (!equiv)  scales[wname] = xsec->point(0).x() / sumw;
        }
        else if (!equiv) {
          throw UserError("Missing cross-section, needed for non-equivalent merging!");
        }
      }
      // Now check if any options should be removed
      for ( const string& delopt : delopts ) {
        if ( path.hasOption(delopt) )  path.removeOption(delopt);
      }
      // ...or added
      for (size_t i = 0; i < optAnas.size(); ++i) {
        if (path.path().find(optAnas[i]) != string::npos ) {
          path.setOption(optKeys[i], optVals[i]);
          path.fixOptionString();
        }
      }
      path.setPath();
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      // merge AOs
      const string& key = path.path();
      std::cout << __FILE__ << ": " << __LINE__ << ": " << key << std::endl;
      const double sf = key.find("_EVTCOUNT") != string::npos? 1 : scales[wname];
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      std::cout << sf << std::endl;
      auto it = allaos.find(key);
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      std::cout << (it == allaos.end()) << std::endl;
      if (allaos.find(key) == allaos.end()) {
        std::cout << __FILE__ << ": " << __LINE__ << std::endl;
        MSG_DEBUG("Copy first occurrence of " << key << " using scale " << sf);
        allaos[key] = ao; // TODO would be nice to combine these two?
        std::cout << __FILE__ << ": " << __LINE__ << std::endl;
        copyao(ao, allaos[key], sf);
        std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      }
      else if ( !addaos(allaos[key], ao, sf) ) {
        std::cout << __FILE__ << ": " << __LINE__ << std::endl;
        MSG_DEBUG("Cannot merge objects with path " << key
                  << " of type " << ao->annotation("Type") << " using scale " << sf);
      } // end of merge attempt
    } // loop over all AOs ends
  }

  void AnalysisHandler::mergeAOS(map<string, YODA::AnalysisObjectPtr> &allaos,
                                 map<string, YODA::AnalysisObjectPtr> &newaos, 
                                 map<string, pair<double, double>> &allxsecs,
                                 const vector<string> &delopts,
                                 const vector<string> &optAnas,
                                 const vector<string> &optKeys,
                                 const vector<string> &optVals,
                                 const bool equiv, 
                                 const bool overwrite_xsec,
                                 const double fileweight) {


    for (const auto & ao : newaos){
      std::cout << ao.first << ", " << ao.second << " (" << ao.second->path() << ")" << std::endl;
    }
                                   
    std::cout << __FILE__ << ": " << __LINE__ << std::endl;
    map<string, double> scales;
    for (auto& [aopath, aor] : newaos) { 
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      YODA::AnalysisObjectPtr ao(aor);
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      //AOPath path(ao->path());
      AOPath path(aopath);
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      if ( !path ) {
        throw UserError("Invalid path name in new AO set!");
      }
      // skip everything that isn't pre-finalize
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      if ( !path.isRaw() ){
        std::cout << __FILE__ << ": " << __LINE__ << std::endl;
        continue;
      } 

      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      MSG_DEBUG(" " << ao->path());

      const string& wname = path.weightComponent();
      if ( scales.find(wname) == scales.end() ) {
        std::cout << __FILE__ << ": " << __LINE__ << std::endl;
        scales[wname] = 1.0;
        // get the sum of weights and number of entries for the current weight
        double evts = 0, sumw = 1;
        auto ec_it = newaos.find("/RAW/_EVTCOUNT" + wname);
        if ( ec_it != newaos.end() ) {
          YODA::Counter* cPtr = static_cast<YODA::Counter*>(ec_it->second.get());
          evts = cPtr->numEntries();
          sumw = cPtr->sumW()? cPtr->sumW() : 1;
        }
        else if (!equiv) {
          throw UserError("Missing event counter, needed for non-equivalent merging!");
        }
        std::cout << __FILE__ << ": " << __LINE__ << std::endl;
        // in stacking mode: add up all the cross sections
        // in equivalent mode: weight the cross-sections
        // estimates by the corresponding number of entries
        const string xspath = "/RAW/_XSEC" + wname;
        auto xs_it = newaos.find(xspath);
        if ( xs_it != newaos.end() ) {
          YODA::Scatter1D* xsec = static_cast<YODA::Scatter1D*>(xs_it->second.get());
          if (overwrite_xsec) {
            MSG_DEBUG("Set user-supplied weight: " << fileweight);
            xsec->point(0).setX(fileweight);
          }
          else {
            MSG_DEBUG("Multiply user-supplied weight: " << fileweight);
            xsec->scaleX(fileweight);
          }
          // get iterator to the existing (or newly created) key-value pair
          auto xit = allxsecs.insert( make_pair(xspath, make_pair(0,0)) ).first;
          // update cross-sections, possibly weighted by number of entries
          xit->second.first  += (equiv? evts : 1.0) * xsec->point(0).x();
          xit->second.second += (equiv? sqr(evts) : 1.0) * sqr(xsec->point(0).xErrAvg());
          // only in stacking mode: multiply each AO by cross-section / sumW
          if (!equiv)  scales[wname] = xsec->point(0).x() / sumw;
        }
        else if (!equiv) {
          throw UserError("Missing cross-section, needed for non-equivalent merging!");
        }
      }
      // Now check if any options should be removed
      for ( const string& delopt : delopts ) {
        if ( path.hasOption(delopt) )  path.removeOption(delopt);
      }
      // ...or added
      for (size_t i = 0; i < optAnas.size(); ++i) {
        if (path.path().find(optAnas[i]) != string::npos ) {
          path.setOption(optKeys[i], optVals[i]);
          path.fixOptionString();
        }
      }
      path.setPath();
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      // merge AOs
      const string& key = path.path();
      std::cout << __FILE__ << ": " << __LINE__ << ": " << key << std::endl;
      const double sf = key.find("_EVTCOUNT") != string::npos? 1 : scales[wname];
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      std::cout << sf << std::endl;
      auto it = allaos.find(key);
      std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      std::cout << (it == allaos.end()) << std::endl;
      if (allaos.find(key) == allaos.end()) {
        std::cout << __FILE__ << ": " << __LINE__ << std::endl;
        MSG_DEBUG("Copy first occurrence of " << key << " using scale " << sf);
        allaos[key] = ao; // TODO would be nice to combine these two?
        std::cout << __FILE__ << ": " << __LINE__ << std::endl;
        copyao(ao, allaos[key], sf);
        std::cout << __FILE__ << ": " << __LINE__ << std::endl;
      }
      else if ( !addaos(allaos[key], ao, sf) ) {
        std::cout << __FILE__ << ": " << __LINE__ << std::endl;
        MSG_DEBUG("Cannot merge objects with path " << key
                  << " of type " << ao->annotation("Type") << " using scale " << sf);
      } // end of merge attempt
    } // loop over all AOs ends
  }


  void AnalysisHandler::setupReentrantRun(map<string, YODA::AnalysisObjectPtr> reentrantAOs,
                                          map<string, pair<double,double> > reentrantXsecs,
                                          const bool equiv) {

    // get list of analyses & multi-weights to be initialised
    set<string> foundAnalyses;
    set<string> foundWeightNames;
    for (const auto& pair : reentrantAOs) {
      AOPath path(pair.first); 
      if ( path.analysisWithOptions() != "" ) {
        foundAnalyses.insert(path.analysisWithOptions());
      }
      foundWeightNames.insert(path.weight());
    }

    // Make analysis handler aware of the weight names present
    _weightNames.clear();
    _rivetDefaultWeightIdx = _defaultWeightIdx = 0;
    _weightNames = vector<string>(foundWeightNames.begin(), foundWeightNames.end());

    // Then we create and initialize all analyses
    for (const string& ananame : foundAnalyses) { addAnalysis(ananame); }
    _stage = Stage::INIT;
    for (AnaHandle a : analyses() ) {
      MSG_TRACE("Initialising analysis: " << a->name());
      if ( !a->info().reentrant() )
        MSG_WARNING("Analysis " << a->name() << " has not been validated to have "
                    << "a reentrant finalize method. The merged result is unpredictable.");
      try {
        // Allow projection registration in the init phase onwards
        a->_allowProjReg = true;
        a->setProjectionHandler(_projHandler);
        a->init();
        a->syncDeclQueue();
      } catch (const Error& err) {
        cerr << "Error in " << a->name() << "::init method: " << err.what() << endl;
        exit(1);
      }
      MSG_TRACE("Done initialising analysis: " << a->name());
    } // analyses
    _stage = Stage::OTHER;
    _initialised = true;

    // Collect global weights and cross sections and fix scaling for all files
    MSG_DEBUG("Getting event counter and cross-section from "
              << weightNames().size() << " " << numWeights());
    _eventCounter = CounterPtr(weightNames(), Counter("_EVTCOUNT"));
    _xs = Scatter1DPtr(weightNames(), Scatter1D("_XSEC"));
    vector<double> scales(numWeights(), 1.0);
    for (size_t iW = 0; iW < numWeights(); ++iW) {
      MSG_DEBUG("Weight # " << iW << " of " << numWeights());
      _eventCounter.get()->setActiveWeightIdx(iW);
      _xs.get()->setActiveWeightIdx(iW);
      YODA::Scatter1D & xsec = *_xs;
      // set the sum of weights
      auto aoit = reentrantAOs.find(_eventCounter->path());
      if (aoit != reentrantAOs.end()) {
        *_eventCounter += *dynamic_pointer_cast<YODA::Counter>(aoit->second);
      }

      const auto xit = reentrantXsecs.find(xsec.path());
      if ( xit != reentrantXsecs.end() ) {
        double xs = xit->second.first;
        double xserr = sqrt(xit->second.second);
        if ( equiv ) {
          MSG_DEBUG("Equivalent mode: scale by numEntries");
          const double nentries = _eventCounter->numEntries();
          xs /= nentries;
          xserr /= nentries;
        }
        else if (xs) {
          // in stacking mode: need to unscale prior to finalize
          scales[iW] = _eventCounter->sumW()/xs;
        }
        xsec.reset();
        xsec.addPoint( Point1D(xs,xserr) );
      }
      else {
        throw UserError("Missing cross-section for " + xsec.path());
      }

      // Go through all analyses and add stuff to their analysis objects;
      for (AnaHandle a : analyses()) {
        for (const auto& ao : a->analysisObjects()) {
          ao.get()->setActiveWeightIdx(iW);
          YODA::AnalysisObjectPtr yao = ao.get()->activeYODAPtr();
          auto aoit = reentrantAOs.find(yao->path());
          if (aoit != reentrantAOs.end()) {
            if ( !addaos(yao, aoit->second, scales[iW]) ) {
              MSG_DEBUG("Overwriting incompatible starting version of " << yao->path()
                        << " using scale " << scales[iW]);
              copyao(aoit->second, yao, 1.0); // input already scaled by addaos
            }
          }
          else {
            MSG_DEBUG("Cannot merge objects with path " << yao->path()
                      << " of type " << yao->annotation("Type"));
          }
          a->rawHookIn(yao);
          ao.get()->unsetActiveWeight();
        }
      }
      _eventCounter.get()->unsetActiveWeight();
      _xs.get()->unsetActiveWeight();
    }
  }


  void AnalysisHandler::readData(const string& filename) {
    try {
      /// @todo Use new YODA SFINAE to fill the smart ptr vector directly
      vector<YODA::AnalysisObject*> aos_raw;
      YODA::read(filename, aos_raw);
      for (YODA::AnalysisObject* aor : aos_raw)
        _preloads[aor->path()] = YODA::AnalysisObjectPtr(aor);
    } catch (...) { //< YODA::ReadError&
      throw UserError("Unexpected error in reading file: " + filename);
    }
  }


  vector<MultiweightAOPtr> AnalysisHandler::getRivetAOs() const {
      vector<MultiweightAOPtr> rtn;

      for (AnaHandle a : analyses()) {
          for (const auto & ao : a->analysisObjects()) {
              rtn.push_back(ao);
          }
      }
      rtn.push_back(_eventCounter);
      rtn.push_back(_xs);
      return rtn;
  }


  vector<YODA::AnalysisObjectPtr> AnalysisHandler::getYodaAOs(bool includeraw) const {
    vector<YODA::AnalysisObjectPtr> output;

    // First get all multiweight AOs
    vector<MultiweightAOPtr> raos = getRivetAOs();
    output.reserve(raos.size() * numWeights() * (includeraw ? 2 : 1));

    // Identify an index ordering so that default weight is written out first
    vector<size_t> order = { _rivetDefaultWeightIdx };
    for ( size_t  i = 0; i < numWeights(); ++i ) {
      if ( i != _rivetDefaultWeightIdx )  order.push_back(i);
    }

    // Then we go through all finalized AOs one weight at a time
    for (size_t iW : order ) {
      for ( auto rao : raos ) {
        rao.get()->setActiveFinalWeightIdx(iW);
        if ( rao->path().find("/TMP/") != string::npos ) continue;
        output.push_back(rao.get()->activeYODAPtr());
      }
    }

    // Analyses can make changes neccesary for merging to RAW objects
    // before writing.
    for (size_t iW : order)
      for (auto a : analyses()) a->rawHookOut(raos, iW);

    // Finally the RAW objects.
    if (includeraw) {
      for (size_t iW : order ) {
        for ( auto rao : raos ) {
          rao.get()->setActiveWeightIdx(iW);
          output.push_back(rao.get()->activeYODAPtr());
        }
      }
    }

    return output;
  }


  void AnalysisHandler::writeData(std::ostream& ostr, const string& fmt) const {

    const vector<YODA::AnalysisObjectPtr> output = getYodaAOs(true);
    try {
      YODA::write(ostr, begin(output), end(output), fmt);
    } catch (...) { //< YODA::WriteError&
      throw UserError("Unexpected error in writing output");
    }

  }


  void AnalysisHandler::writeData(const string& filename) const {

    const vector<YODA::AnalysisObjectPtr> output = getYodaAOs(true);
    try {
      YODA::write(filename, begin(output), end(output));
    } catch (...) { //< YODA::WriteError&
      throw UserError("Unexpected error in writing file: " + filename);
    }

  }


  string AnalysisHandler::runName() const {
    return _runname;
  }


  size_t AnalysisHandler::numEvents() const {
    return _eventCounter->numEntries();
  }


  std::vector<std::string> AnalysisHandler::analysisNames() const {
    std::vector<std::string> rtn;
    for (AnaHandle a : analyses()) {
      rtn.push_back(a->name());
    }
    return rtn;
  }


  std::vector<std::string> AnalysisHandler::stdAnalysisNames() const {
    // std::vector<std::string> rtn;
    // const string anadatpath = findAnalysisDataFile("analyses.dat");
    // if (fileexists(anadatpath)) {
    //   std::ifstream anadat(anadatpath);
    //   string ananame;
    //   while (anadat >> ananame) rtn += ananame;
    // }
    // return rtn;
    return AnalysisLoader::stdAnalysisNames();
  }


  AnalysisHandler& AnalysisHandler::addAnalyses(const std::vector<std::string>& analysisnames) {
    for (const string& aname : analysisnames) {
      //MSG_DEBUG("Adding analysis '" << aname << "'");
      addAnalysis(aname);
    }
    return *this;
  }


  AnalysisHandler& AnalysisHandler::removeAnalyses(const std::vector<std::string>& analysisnames) {
    for (const string& aname : analysisnames) removeAnalysis(aname);
    return *this;
  }


  void AnalysisHandler::setCrossSection(const pair<double,double>& xsec, bool isUserSupplied) {
    // Update the user xsec
    if (isUserSupplied) {
      MSG_DEBUG("Setting user cross-section = " << xsec.first << " +- " << xsec.second << " pb");
      _userxs = xsec;
    }

    // If not setting the user xsec, and a user xsec is already set, do nothing and exit early
    if (!isUserSupplied && notNaN(_userxs.first)) return;

    // Otherwise, update the xs scatter: xs_var = xs_nom * (sumW_var/sumW_nom)
    /// @todo Performance optimization? Overwriting the whole scatter wrapper on every event seems inefficient...
    MSG_TRACE("Setting nominal cross-section = " << xsec.first << " +- " << xsec.second << " pb");
    _xs = Scatter1DPtr(weightNames(), Scatter1D("_XSEC"));
    _eventCounter.get()->setActiveWeightIdx(_rivetDefaultWeightIdx);
    const double nomwgt = sumW();
    const double nomwt2 = sumW2();
    for (size_t iW = 0; iW < numWeights(); ++iW) {
      _eventCounter.get()->setActiveWeightIdx(iW);
      const double s  = nomwgt? (sumW() / nomwgt) : 1.0;
      const double s2 = nomwt2? sqrt(sumW2() / nomwt2) : 1.0;
      _xs.get()->setActiveWeightIdx(iW);
      _xs->addPoint(xsec.first*s, xsec.second*s2);
    }
    _eventCounter.get()->unsetActiveWeight();
    _xs.get()->unsetActiveWeight();
  }


  double AnalysisHandler::nominalCrossSection() const {
    _xs.get()->setActiveWeightIdx(_rivetDefaultWeightIdx);
    const YODA::Scatter1D::Points& ps = _xs->points();
    if (ps.size() != 1) {
      string errMsg = "Value missing when requesting nominal cross-section";
      throw Error(errMsg);
    }
    double xs = ps[0].x();
    _xs.get()->unsetActiveWeight();
    return xs;
  }


  AnalysisHandler& AnalysisHandler::addAnalysis(Analysis* analysis) {
    analysis->_analysishandler = this;
    // _analyses.insert(AnaHandle(analysis));
    _analyses[analysis->name()] = AnaHandle(analysis);
    return *this;
  }



  AnalysisHandler& AnalysisHandler::setRunBeams(const ParticlePair& beams) {
    _beams = beams;
    MSG_DEBUG("Setting run beams = " << beams << " @ " << sqrtS(beams)/GeV << " GeV");
    return *this;
  }


  PdgIdPair AnalysisHandler::runBeamIDs() const {
    return pids(runBeams());
  }

  pair<double,double> AnalysisHandler::runBeamEnergies() const {
    return energies(runBeams());
  }

  double AnalysisHandler::runSqrtS() const {
    return sqrtS(runBeams());
  }


//TODO @TP: It would be nifty and more efficient if we could merge multiple handlers at once using variadic templates (come back to when it works)
// Merges AnalysisHandler Other into this.
// Not (yet?) written to be symmetric: looks to merge from other into this: if other has stuff (e.g. an analysis)
// that this does not have, it will be ignored. VERY VERY VERY much a WIP.
//TODO: Include weights here?
  void AnalysisHandler::mergeAnalysisHandlers(AnalysisHandler& other, bool equiv){
    MSG_TRACE("Merging analysis handler " << &other << " into " << this);

    //Handlers to be merged must have same beam:
    //TODO: Would it make sense to have a Rivet particle operator== 
    if (other._beams.first.energy() != _beams.first.energy() &&
         other._beams.second.energy() != _beams.second.energy()){return;}

    //TODO @TP: Should we check the "finalisation status" of analysishandlers? If one is and one isn't,
    //things could go weird - is the "stage" enum "mature" enough for this use?

    //Sum event numbers and the YODA counter
    _eventNumber = _eventNumber + other._eventNumber;
    //TODO @TP: This looks horrible, is there a less stomach-churning syntax?
    MSG_TRACE("This evtcounter: " << _eventCounter->val());
    MSG_TRACE("Other evtcounter: " << other._eventCounter->val());
    _eventCounter->operator+=(*(other._eventCounter));
    MSG_TRACE("After merging, this evtcounter: " << _eventCounter->val());

    //TODO @TP: what to do about XS? (I've messed this up like twice now)
    //MSG_TRACE("This XS: " << _xs->


    //Let's look at the weights.
    std::string weightstring; for (auto w : _weightNames){weightstring+=std::string(w+", ");}
    std::string otherweightstring; for (auto w : other._weightNames){otherweightstring+=std::string(w+", ");}
    MSG_TRACE("This ah's weights (size "<<_weightNames.size() <<"):" << weightstring);
    MSG_TRACE("Other ah's weights (size "<<other._weightNames.size() <<"):" << otherweightstring);
    std::string weightindstring; for (auto w : _weightIndices){weightindstring+=(std::to_string(w)+", ");}
    std::string otherweightindstring; for (auto w : other._weightIndices){otherweightindstring+=(std::to_string(w)+", ");}
    MSG_TRACE("This ah's weight indices (size"<<_weightIndices.size()<< "):" << weightindstring);
    MSG_TRACE("Other ah's weight indices (size"<<other._weightIndices.size()<<"):" << otherweightindstring);

    MSG_TRACE("This ah's (_matchWN, _unmatchWN, _nominalWeightNames): ("<<_matchWeightNames<< ", "<<_unmatchWeightNames<<", "<<_nominalWeightName<<")");
    MSG_TRACE("Other ah's (_matchWN, _unmatchWN, _nominalWeightNames): ("<<other._matchWeightNames<< ", "<<other._unmatchWeightNames<<", "<<other._nominalWeightName<<")");

    MSG_TRACE("This ah's subEventWeightsSize is " << _subEventWeights.size());
    MSG_TRACE("Other ah's subEventWeightsSize is " << other._subEventWeights.size());
    std::cout << "SubEventWeights: (";
    for(auto i : _subEventWeights){
      std::cout << "(";
      for (auto j : i){
        std::cout << j << ", ";
      }
      std::cout << "), ";
    }
    std::cout << ")\nOther subeventweights: (";
    for(auto i : other._subEventWeights){
      std::cout << "(";
      for (auto j : i){
        std::cout << j << ", ";
      }
      std::cout << "), ";
    }
    std::cout << ")\n";



    size_t iW = 1;

    const std::vector<AnaHandle> othersAnalyses = other.analyses();

    for(AnaHandle a : analyses()){
      //find corresponding analysis in the other ao;
      std::string analysisname = a->name();
      auto other_analysis_it = std::find_if(othersAnalyses.begin(), othersAnalyses.end(),
                   [&analysisname](const AnaHandle otherhandle){std::cerr<<"\n"<<otherhandle->name() <<" v "<< analysisname<<endl;return otherhandle->name() == analysisname; });

      if (other_analysis_it == othersAnalyses.end()){
        //This analysis isn't present in the other ah. Move on.
        MSG_DEBUG("Analysis " << a->name() << " present in analysishandler " << this 
                      << " not present in " << &other);
        //Debug only:
        // std::string analysesinother;
        // for (auto ana : other.analyses()){
        //   analysesinother += std::string(ana->name() + ", ");
        // }
        // MSG_TRACE("Analyses present in other: " << analysesinother);
        continue;
      }

      //TODO: Is there a more elegant syntax for double dereferencing? This looks awful.
      auto othersaos = (*other_analysis_it)->analysisObjects();
      //Is it a valid assumption that the analysisObjects() only contains one object of each name? I assume so?
      MSG_TRACE("There are " << othersaos.size() << " other aos");
      std::string otheraostring; for (auto ao : othersaos){otheraostring+=std::string(ao->name()+", ");}
      MSG_TRACE("They are " << otheraostring);
      MSG_TRACE("This analysis has " << a->analysisObjects().size() << " aos");
      std::string aostring; for (auto ao : a->analysisObjects()){aostring+=std::string(ao->name()+", ");}
      MSG_TRACE("They are " << aostring);
      for (const auto& ao : a-> analysisObjects()){
        //Find corresponding ao in other:
        auto other_ao_it = std::find_if(othersaos.begin(), othersaos.end(),
                                        [&](const Rivet::MultiweightAOPtr otherAOptr){return ao->name() == otherAOptr->name();});
        if (other_ao_it == (*other_analysis_it)->analysisObjects().end()){
          //This analysis object is not present in other. Slightly scary. Move on (but warn?!)
          MSG_WARNING("Analysis object " << ao->name() << " present in analysis " << a->name() << 
                "  of analysishandler " << this << " NOT found in same analysis of analysishandler "
                << &other << ". AnalysisHandler merge has probably gone wrong..." );
          continue;
        }
        MSG_TRACE("Merging ao " << ao->name() << " in analysis " << a->name());

        //ao.get()->setActiveWeightIdx(iW);
        YODA::AnalysisObjectPtr yao = ao.get()->activeYODAPtr();
        if (!addaos(yao, other_ao_it->get()->activeYODAPtr(), iW )){
          MSG_WARNING("Failed to merge object named " << ao->name() << " in analysis " << a->name()
                      << " for analysis handlers " << this << " and " << &other << 
                      " (if this ao is not a histogram, this is nothing to worry about)");
        }
        //ao.get()->unsetActiveWeight();

      }

    }

    MSG_TRACE("DEBUG INFO");
    print_ah_info();


    //Finalise all analyses;
    //TODO -> do I actually want this here?
    MSG_TRACE("Before calling finalize, this evtcounter: " << _eventCounter->val());
    finalize();
    MSG_TRACE("After calling finalize, this evtcounter: " << _eventCounter->val());
  }

  //TODO: I think there's potential to vectorise this.
  void AnalysisHandler::combineAnalysisHandlers(std::vector<AnalysisHandler*> &handlers, bool equiv){
    //Variables that will go into mergeAOs
    std::map<string, YODA::AnalysisObjectPtr> allaos;
    std::map<string, YODA::AnalysisObjectPtr> newaos;
    map<string, pair<double, double>> allxsecs;

    std::cout << __FILE__ << ": " << __LINE__ << std::endl;

    for (const AnalysisHandler* handler: handlers){
      for (const YODA::AnalysisObjectPtr& yodaptr : handler->getYodaAOs(true)){
        AOPath aopath(yodaptr->path());
        if (aopath.isRaw()){
          newaos[yodaptr->path()] = yodaptr;
        }
      }
      mergeAOS(allaos, newaos, allxsecs, {}, {}, {}, {}, equiv);
    }

    std::cout << __FILE__ << ": " << __LINE__ << std::endl;
    setupReentrantRun(allaos, allxsecs, equiv);
    return;
  }

}
