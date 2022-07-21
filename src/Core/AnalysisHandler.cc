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

using std::cout;
using std::cerr;

namespace Rivet {


  AnalysisHandler::AnalysisHandler(const string& runname)
    : _runname(runname), _userxs{NAN, NAN},
      _initialised(false), _ignoreBeams(false),
      _skipWeights(false), _matchWeightNames(""),
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

    /// @todo Should the Rivet analysis objects know about weight names?

    /// @todo Get the HepMC3::GenRunInfo object from the first event and store/wrap it?

    setRunBeams(Rivet::beams(ge));
    MSG_DEBUG("Initialising the analysis handler");
    _eventNumber = ge.event_number();

    // Assemble the weight streams to be used
    setWeightNames(ge);
    if (_skipWeights) {
        MSG_INFO("Only using nominal weight. Variation weights will be ignored.");
    } else if (haveNamedWeights()) {
        MSG_INFO("Using named weights");
    } else {
        MSG_INFO("NOT using named weights. Using first weight as nominal weight");
    }

    // Create the multi-weighted event counter
    _eventCounter = CounterPtr(weightNames(), Counter("_EVTCOUNT"));

    // Set the cross section based on what is reported by the init-event, else zero
    if (ge.cross_section()) {
      setCrossSection(HepMCUtils::crossSection(ge, _defaultWeightIdx));
    } else {
      MSG_DEBUG("No cross-section detected in first event: setting default to 0 pb");
      setCrossSection({0.0, 0.0});
    }

    // Check that analyses are beam-compatible, and note those that aren't
    const size_t num_anas_requested = analysisNames().size();
    vector<string> anamestodelete;
    for (const AnaHandle& a : analyses()) {
      if (!_ignoreBeams && !a->isCompatible(beams())) {
        //MSG_DEBUG(a->name() << " requires beams " << a->requiredBeams() << " @ " << a->requiredEnergies() << " GeV");
        anamestodelete.push_back(a->name());
      }
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
        a->init();
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
    string nom_winner = "";
    vector<string> nom_shortlist;
    _weightIndices.clear();
    for (size_t i = 0, N = _weightNames.size(); i < N; ++i) {
      _weightIndices.push_back(i);
      if (_weightNames[i] == _nominalWeightName) {
        nom_shortlist.push_back("'"+ _weightNames[i] +"'");
        if (nDefaults == 0) {
          nom_winner = "'" + _weightNames[i] + "'";
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
          nom_shortlist.push_back("'"+ _weightNames[i] +"'");
          if (nDefaults == 0) {
            nom_winner = "'" + _weightNames[i] + "'";
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
      MSG_WARNING("Found " << nDefaults << " default weight candidates: " << join(nom_shortlist, ", ") << ". Will use: " <<  nom_winner);
    }

    // Apply behaviours for only using the nominal weight, or all weights
    if (_skipWeights)  {

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
            MSG_DEBUG("Selected nominal weight: " << _weightNames[i]);
            continue;
          }
          for (const std::regex& re : patterns) {
            if ( std::regex_match(_weightNames[i], re) ) {
              selected_indices.push_back(_weightIndices[i]);
              selected_subset.push_back(_weightNames[i]);
              MSG_DEBUG("Selected variation weight: " << _weightNames[i]);
              break;
            }
          }
        }
        _weightNames = selected_subset;
        _weightIndices = selected_indices;
      }

      // Check if the remaining weight names match supplied string/regexes and *de*select accordingly
      vector<std::regex> patterns = { std::regex("^IRREG.*", std::regex_constants::icase) };
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


  void AnalysisHandler::analyze(const GenEvent& ge) {
    // Call init with event as template if not already initialised
    if (!_initialised) init(ge);
    assert(_initialised);

    // Ensure that beam details match those from the first event (if we're checking beams)
    if ( !_ignoreBeams ) {
      const PdgIdPair evtbeams = Rivet::beamIds(ge);
      const double sqrts = Rivet::sqrtS(ge);
      MSG_DEBUG("Event beams = " << evtbeams << " at sqrt(s) = " << sqrts/GeV << " GeV");
      if (evtbeams.first == PID::ANY && evtbeams.second == PID::ANY) {
        MSG_ERROR("No event beams found: please fix the events, or run with beam-checking disabled");
        exit(1);
      }
      if (!compatible(evtbeams, _beams) || !fuzzyEquals(sqrts, sqrtS())) {
        cerr << "Event beams mismatch: "
             << PID::toBeamsString(evtbeams) << " @ " << sqrts/GeV << " GeV" << " vs. first beams "
             << this->beams() << " @ " << this->sqrtS()/GeV << " GeV" << endl;
        exit(1);
      }
    }

    // Create the Rivet event wrapper
    /// @todo Filter/normalize the event here
    /// @todo Find a way to cache the env call
    bool strip = ( getEnvParam("RIVET_STRIP_HEPMC", string("NOOOO") ) != "NOOOO" );
    Event event(ge, _weightIndices, strip);

    // Set the cross section based on what is reported by this event
    if (ge.cross_section())  setCrossSection(event.crossSections());

    // If the event number has changed, sync the sub-event analysis objects to persistent
    // NB. Won't happen for first event because _eventNumber is set in init()
    if (_eventNumber != ge.event_number()) {
      pushToPersistent();
      _eventNumber = ge.event_number();

      // Dump current final histograms
      if ( _dumpPeriod > 0 && numEvents() > 0 && numEvents() % _dumpPeriod == 0 ) {
        MSG_DEBUG("Dumping intermediate results to " << _dumpFile << ".");
        _dumping = numEvents()/_dumpPeriod;
        finalize();
        writeData(_dumpFile);
        _dumping = 0;
      }

    }

    // Make a new sub-event: affects every analysis object
    MSG_TRACE("Starting new sub-event");
    _eventCounter.get()->newSubEvent();
    for (const AnaHandle& a : analyses()) {
      for (auto ao : a->analysisObjects()) {
        ao.get()->newSubEvent();
      }
    }

    // Append to sub-event weights list, modulo weight-capping to avoid spikes
    _subEventWeights.push_back(event.weights());
    if (_weightCap != 0.) {
      MSG_DEBUG("Implementing weight cap using a maximum |weight| = " << _weightCap << " for latest subevent.");
      size_t lastSub = _subEventWeights.size() - 1;
      for (size_t i = 0; i < _subEventWeights[lastSub].size(); ++i) {
        if (abs(_subEventWeights[lastSub][i]) > _weightCap) {
          _subEventWeights[lastSub][i] = sign(_subEventWeights[lastSub][i]) * _weightCap;
        }
      }
    }
    MSG_DEBUG("Analyzing subevent #" << _subEventWeights.size() - 1 << ".");

    // Warn if the subevent list is getting very long without flushing
    if (_subEventWeights.size() % 1000 == 0) {
      MSG_WARNING("Sub-event weight list has " << _subEventWeights.size() << " elements: are the weight numbers correctly set in the input events?");
    }

    // Update the event counter
    // NB. updated on sub-events, but synced cf. histos on full-event boundaries
    _eventCounter->fill();

    // Run the analyses
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

  }


  void AnalysisHandler::analyze(const GenEvent* ge) {
    if (ge == nullptr) {
      MSG_ERROR("AnalysisHandler received null pointer to GenEvent");
      //throw Error("AnalysisHandler received null pointer to GenEvent");
    }
    analyze(*ge);
  }


  void AnalysisHandler::pushToPersistent() {
    if ( _subEventWeights.empty() ) return;
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

    // First push all analyses' objects to persistent and final
    MSG_TRACE("AnalysisHandler::finalize(): Pushing analysis objects to persistent.");
    pushToPersistent();

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
    _eventCounter.get()->setActiveFinalWeightIdx(defaultWeightIndex());
    _xs.get()->setActiveFinalWeightIdx(defaultWeightIndex());
    if (!_dumping) {
      const int nevts = numEvents();
      MSG_DEBUG("Processed " << nevts << " event" << (nevts != 1 ? "s" : ""));
    }

    _stage = Stage::OTHER;

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


  void AnalysisHandler::mergeYodas(const vector<string> &aofiles,
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
        MSG_WARNING("Malformed analysis option: "+addopt+". Format as :OPT=VAL or ANA:OPT=VAL");
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
          // skip everything that isn't pre-finalize
          if ( !AOPath(aopath).isRaw() ) continue;
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
      if (equiv && 2*raw_map.size() != aos_raw.size()) {
        MSG_WARNING("Number of pre- and post-finalize AOs do not match for file: " << file);
        continue;
      }

      // merge AOs from current file into "allaos"
      mergeAOS(allaos, raw_map, allxsecs, delopts, optAnas, optKeys, optVals,
                                equiv, overwrite_xsec, fileweight);

    } // loop over all input files ends
    std::cout << std::endl;

    if (allaos.empty()) {
      cerr << "Insuffient pre-finalize AOs to do a reentrant run!" << endl;
      exit(1);
    }

    MSG_DEBUG("Finalize cross-section scaling ...");
    vector<double> scales(allxsecs.size(), 1.0);
    for (const auto& item : allxsecs) {
      const string wname = item.first;
      double xs = item.second.first;
      double xserr = sqrt(item.second.second);
      auto xs_it = allaos.find("/RAW/_XSEC" + wname);
      assert( xs_it != allaos.end() );
      shared_ptr<YODA::Scatter1D> xsec = dynamic_pointer_cast<YODA::Scatter1D>(xs_it->second);
      auto ec_it = allaos.find("/RAW/_EVTCOUNT" + wname);
      assert( ec_it != allaos.end() );
      if (equiv) {
        MSG_DEBUG("Equivalent mode: scale by numEntries");
        const double nentries = dynamic_pointer_cast<YODA::Counter>(ec_it->second)->numEntries();
        xs /= nentries;
        xserr /= nentries;
      }
      xsec->reset();
      xsec->addPoint( Point1D(xs,xserr) );
    }

    MSG_INFO("Rerunning finalize ...");

    // initialise analyses and load merged AOs back into memory
    // set unscale to true if equiv is false
    loadAOs(allaos, !equiv);

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
                                 const double user_xsec) {


    map<string, double> scales;
    for (const auto& item : newaos) {
      const string& aopath = item.first;
      YODA::AnalysisObjectPtr ao(item.second);
      //AOPath path(ao->path());
      AOPath path(aopath);
      if ( !path ) {
        throw UserError("Invalid path name in new AO set!");
      }
      // skip everything that isn't pre-finalize
      if ( !path.isRaw() ) continue;

      MSG_DEBUG(" " << ao->path());

      const string& wname = path.weightComponent();
      if ( scales.find(wname) == scales.end() ) {
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
        // in stacking mode: add up all the cross sections
        // in equivalent mode: weight the cross-sections
        // estimates by the corresponding number of entries
        const string xspath = "/RAW/_XSEC" + wname;
        auto xs_it = newaos.find(xspath);
        if ( xs_it != newaos.end() ) {
          YODA::Scatter1D* xsec = static_cast<YODA::Scatter1D*>(xs_it->second);
          if (overwrite_xsec) {
            MSG_DEBUG("Set user-supplied weight: " << user_xsec);
            xsec->point(0).setX(user_xsec);
          }
          else {
            MSG_DEBUG("Multiply user-supplied weight: " << user_xsec);
            xsec->scaleX(user_xsec);
          }
          // get iterator to the existing (or newly created) key-value pair
          auto xit = allxsecs.insert( make_pair(wname, make_pair(0,0)) ).first;
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
        if (optAnas[i] == "" || path.path().find(optAnas[i]) != string::npos ) {
          path.setOption(optKeys[i], optVals[i]);
          path.fixOptionString();
        }
      }
      path.setPath();

      // merge AOs
      const string& key = path.path();
      const double sf = key.find("_EVTCOUNT") != string::npos? 1 : scales[wname];
      if (allaos.find(key) == allaos.end()) {
        MSG_DEBUG("Copy first occurrence of " << key << " using scale " << sf);
        allaos[key] = ao; // TODO would be nice to combine these two?
        copyao(ao, allaos[key], sf);
      }
      else if ( !addaos(allaos[key], ao, sf) ) {
        MSG_DEBUG("Cannot merge objects with path " << key
                  << " of type " << ao->annotation("Type") << " using scale " << sf);
      } // end of merge attempt
    } // loop over all new AOs ends
  }


  void AnalysisHandler::loadAOs(const map<string, YODA::AnalysisObjectPtr>& allAOs, const bool unscale) {

    // Check that AH hasn't already been initialised
    if (_initialised)
      throw UserError("AnalysisHandler::init has already been called: cannot re-initialize!");


    // get list of analyses & multi-weights to be initialised
    set<string> foundAnalyses;
    set<string> foundWeightNames;
    for (const auto& pair : allAOs) {
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
        a->init();
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
      // set the sum of weights
      auto aoit = allAOs.find(_eventCounter->path());
      if (aoit != allAOs.end()) {
        *_eventCounter = *dynamic_pointer_cast<YODA::Counter>(aoit->second);
      }

      // set the cross-section
      const auto xit = allAOs.find(_xs->path());
      if ( xit != allAOs.end() ) {
        *_xs = *dynamic_pointer_cast<YODA::Scatter1D>(xit->second);
        if (unscale && _xs->point(0).x()) {
          // in stacking mode: need to unscale prior to finalize
          scales[iW] = _eventCounter->sumW()/_xs->point(0).x();
        }
      }
      else {
        throw UserError("Missing cross-section for " + _xs->path());
      }

      // Go through all analyses and add stuff to their analysis objects;
      for (AnaHandle a : analyses()) {
        for (const auto& ao : a->analysisObjects()) {
          ao.get()->setActiveWeightIdx(iW);
          YODA::AnalysisObjectPtr yao = ao.get()->activeYODAPtr();
          auto aoit = allAOs.find(yao->path());
          if (aoit != allAOs.end()) {
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


  void AnalysisHandler::merge(AnalysisHandler& other) {

    // Check if both AHs have been initialised
    if (!_initialised || !other._initialised)
      throw UserError("AnalysisHandler::init has not been called: cannot merge!");

    // Check if both AHs contain the same registered analyses
    const std::vector<std::string> &this_anaNames = analysisNames();
    const std::vector<std::string> &that_anaNames = other.analysisNames();
    bool is_equal = bool(this_anaNames.size() == that_anaNames.size());
    if (is_equal)   is_equal = std::equal(this_anaNames.begin(), this_anaNames.end(), that_anaNames.begin());
    if (!is_equal)  throw UserError("The AnalysisHandlers are not equivalent!");

    // @todo Do we need to check that the sequence of weight indices is the same?

    // Check if the registered analyses are reentrant safe
    for (AnaHandle a : analyses() ) {
      MSG_TRACE("Initialising analysis: " << a->name());
      if ( !a->info().reentrant() )
        MSG_WARNING("Analysis " << a->name() << " has not been validated to have "
                    << "a reentrant finalize method. The merged result is unpredictable.");
    }
    _stage = Stage::OTHER;

    // First push all analyses' objects to persistent and final
    MSG_TRACE("AnalysisHandler::merge(): Pushing analysis objects to persistent.");
    pushToPersistent();
    other.pushToPersistent();

    // Collect global weights and cross sections and fix scaling for all AHs
    MSG_DEBUG("Getting event counter and cross-section from "
              << weightNames().size() << " " << numWeights());

    for (size_t iW = 0; iW < numWeights(); ++iW) {
      MSG_DEBUG("Weight # " << iW << " of " << numWeights());
      // set the sum of weights
      _eventCounter.get()->setActiveWeightIdx(iW);
      const double this_evts = _eventCounter->numEntries();
      other._eventCounter.get()->setActiveWeightIdx(iW);
      const double other_evts = other._eventCounter->numEntries();
      *_eventCounter += *other._eventCounter;
      // set the cross-section
      _xs.get()->setActiveWeightIdx(iW);
      other._xs.get()->setActiveWeightIdx(iW);
      double xs = this_evts * _xs->point(0).x() + other_evts * other._xs->point(0).x();
      double xserr = sqr(this_evts * _xs->point(0).xErrAvg());
      xserr += sqr(other_evts * other._xs->point(0).xErrAvg());
      const double ntot = _eventCounter->numEntries();
      if (ntot) {
        xs = xs / ntot;
        xserr = sqrt(xserr) / ntot;
      }
      YODA::Scatter1D& this_xs = *_xs;
      this_xs.reset();
      this_xs.addPoint( Point1D(xs,xserr) );


      // Go through all analyses and merge other's AOs with current AOs
      for (const auto& apair : other.analysesMap()) {
        for (const auto& other_ao : apair.second->analysisObjects()) {
          other_ao.get()->setActiveWeightIdx(iW);
          YODA::AnalysisObjectPtr other_yao = other_ao.get()->activeYODAPtr();
          // Find corresponding YODA::AO in current AH
          for (const MultiweightAOPtr& this_ao : analysis(apair.first)->analysisObjects()) {
            this_ao.get()->setActiveWeightIdx(iW);
            if (this_ao->path() != other_ao->path())  continue;
            YODA::AnalysisObjectPtr this_yao = this_ao.get()->activeYODAPtr(); // found it!
            // attempt merge
            if ( !addaos(this_yao, other_yao, 1.0) ) {
              MSG_DEBUG("Overwriting incompatible starting version of " << this_yao->path());
              copyao(other_yao, this_yao, 1.0); // input already scaled by addaos
            }
            analysis(apair.first)->rawHookIn(this_yao);
            this_ao.get()->unsetActiveWeight();
          }
          other_ao.get()->unsetActiveWeight();
          // @todo warn if AO could not be found? Throw an error even?
          //e.g. throw LookupError("Data object " + other_ao->path() + " not found");
        }
      }
      _eventCounter.get()->unsetActiveWeight();
      _xs.get()->unsetActiveWeight();
    } // end of loop over weights
    // leave it to user to call finalize()
  }


  void AnalysisHandler::readData(std::istream& istr, const string& fmt, bool preload) {

    vector<YODA::AnalysisObject*> aos_raw;
    map<string,YODA::AnalysisObjectPtr> aomap;
    try {
      YODA::read(istr, aos_raw, fmt);
      for (YODA::AnalysisObject* aor : aos_raw)
        aomap[aor->path()] = YODA::AnalysisObjectPtr(aor);
    } catch (...) { //< YODA::WriteError&
      throw UserError("Unexpected error in reading input");
    }
    if (preload)  _preloads = std::move(aomap);
    else          loadAOs(aomap);

  }

  void AnalysisHandler::readData(const string& filename, bool preload) {
    map<string,YODA::AnalysisObjectPtr> aomap;
    try {
      /// @todo Use new YODA SFINAE to fill the smart ptr vector directly
      vector<YODA::AnalysisObject*> aos_raw;
      YODA::read(filename, aos_raw);
      for (YODA::AnalysisObject* aor : aos_raw)
        aomap[aor->path()] = YODA::AnalysisObjectPtr(aor);
    } catch (...) { //< YODA::ReadError&
      throw UserError("Unexpected error in reading file: " + filename);
    }
    if (preload)  _preloads = std::move(aomap);
    else          loadAOs(aomap);
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


  void AnalysisHandler::setCrossSection(const vector<pair<double,double>>& xsecs, bool isUserSupplied) {

    if (xsecs.empty())
      throw UserError("No cross-section supplied!");

    if (xsecs.size() == 1)  setCrossSection(xsecs[0], isUserSupplied);
    else {
      // Update the user xsec
      if (isUserSupplied) _userxs = xsecs[0];

      // If not setting the user xsec, and a user xsec is already set, exit early
      if (!isUserSupplied && notNaN(_userxs.first)) return;

      // Otherwise, update the xs scatter
      _xs = Scatter1DPtr(weightNames(), Scatter1D("_XSEC"));
      for (size_t iW = 0; iW < numWeights(); ++iW) {
        _xs.get()->setActiveWeightIdx(iW);
        _xs->addPoint(xsecs[iW].first, xsecs[iW].second);
      }
      _xs.get()->unsetActiveWeight();
    }
  }


  void AnalysisHandler::setCrossSection(const pair<double,double>& xsec, bool isUserSupplied) {
    // Update the user xsec
    if (isUserSupplied) _userxs = xsec;

    // If not setting the user xsec, and a user xsec is already set, exit early
    if (!isUserSupplied && notNaN(_userxs.first)) return;

    // Otherwise, update the xs scatter: xs_var = xs_nom * (sumW_var/sumW_nom)
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


  AnalysisHandler& AnalysisHandler::addAnalysis(Analysis* analysis) {
    analysis->_analysishandler = this;
    // _analyses.insert(AnaHandle(analysis));
    _analyses[analysis->name()] = AnaHandle(analysis);
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


  void AnalysisHandler::skipMultiWeights(bool ignore) {
    _skipWeights = ignore;
  }

  void AnalysisHandler::selectMultiWeights(std::string patterns) {
    _matchWeightNames = patterns;
  }

  void AnalysisHandler::deselectMultiWeights(std::string patterns) {
    _unmatchWeightNames = patterns;
  }

  void AnalysisHandler::setNominalWeightName(std::string name) {
    _nominalWeightName = name;
  }

}
