// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/ParticleName.hh"
#include "Rivet/Tools/BeamConstraint.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Beam.hh"
#include "YODA/ReaderYODA.h"
#include "YODA/WriterYODA.h"
#include <regex>
#include <iostream>

using std::cout;
using std::cerr;

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
      _eventCounter(vector<string>(), Counter()), _xs(vector<string>(), Scatter1D()),
      _initialised(false), _ignoreBeams(false),
      _defaultWeightIdx(0)
  {}


  AnalysisHandler::~AnalysisHandler()
  {  }


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

    /// @todo Should the Rivet analysis objects know about weight names?

    setRunBeams(Rivet::beams(ge));
    MSG_DEBUG("Initialising the analysis handler");
    _eventNumber = ge.event_number();

    setWeightNames(ge);
    if (haveNamedWeights())
        MSG_INFO("Using named weights");
    else
        MSG_INFO("NOT using named weights. Using first weight as nominal weight");

    _eventCounter = CounterPtr(weightNames(), Counter("_EVTCOUNT"));

    // Set the cross section based on what is reported by this event.
    if (ge.cross_section()) {
      MSG_TRACE("Getting cross section.");
      double xs = ge.cross_section()->cross_section();
      double xserr = ge.cross_section()->cross_section_error();
      setCrossSection(xs, xserr);
    }

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
    _stage = Stage::INIT;
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
    _stage = Stage::OTHER;
    _initialised = true;
    MSG_DEBUG("Analysis handler initialised");
  }


  void AnalysisHandler::setWeightNames(const GenEvent& ge) {
    /// reroute the print output to a std::stringstream and process
    /// The iteration is done over a map in hepmc2 so this is safe
    std::ostringstream stream;
    ge.weights().print(stream);  // Super lame, I know
    string str =  stream.str();

    std::regex re("(([^()]+))"); // Regex for stuff enclosed by parentheses ()
    size_t idx = 0;
    for (std::sregex_iterator i = std::sregex_iterator(str.begin(), str.end(), re); i != std::sregex_iterator(); ++i ) {
      std::smatch m = *i;
      vector<string> temp = ::split(m.str(), "[,]");
      if (temp.size() ==2) {
        MSG_DEBUG("Name of weight #" << _weightNames.size() << ": " << temp[0]);

        // store the default weight based on weight names
        if (temp[0] == "Weight" || temp[0] == "0" || temp[0] == "Default") {
          MSG_DEBUG(_weightNames.size() << " is being used as the nominal.");
          _weightNames.push_back("");
          _defaultWeightIdx = idx;
        } else
          _weightNames.push_back(temp[0]);


        idx++;
      }
    }
  }


  void AnalysisHandler::analyze(const GenEvent& ge) {
    // Call init with event as template if not already initialised
    if (!_initialised) init(ge);
    assert(_initialised);

    // Ensure that beam details match those from the first event (if we're checking beams)
    if ( !_ignoreBeams ) {
      const PdgIdPair beams = Rivet::beamIds(ge);
      const double sqrts = Rivet::sqrtS(ge);
      if (!compatible(beams, _beams) || !fuzzyEquals(sqrts, sqrtS())) {
        cerr << "Event beams mismatch: "
             << PID::toBeamsString(beams) << " @ " << sqrts/GeV << " GeV" << " vs. first beams "
             << this->beams() << " @ " << this->sqrtS()/GeV << " GeV" << endl;
        exit(1);
      }
    }


    // Create the Rivet event wrapper
    /// @todo Filter/normalize the event here
    Event event(ge);

    // set the cross section based on what is reported by this event.
    // if no cross section
    MSG_TRACE("getting cross section.");
    if (ge.cross_section()) {
      MSG_TRACE("getting cross section from GenEvent.");
      double xs = ge.cross_section()->cross_section();
      double xserr = ge.cross_section()->cross_section_error();
      setCrossSection(xs, xserr);
    }
    #endif

    // won't happen for first event because _eventNumber is set in
    // init()
    if (_eventNumber != ge.event_number()) {
        /// @todo
        /// can we get away with not passing a matrix?

        MSG_TRACE("AnalysisHandler::analyze(): Pushing _eventCounter to persistent.");
        _eventCounter.get()->pushToPersistent(_subEventWeights);
        // if this is indeed a new event, push the temporary
        // histograms and reset
        for (const AnaHandle& a : _analyses) {
            for (auto ao : a->analysisObjects()) {
                MSG_TRACE("AnalysisHandler::analyze(): Pushing " << a->name() << "'s " << ao->name() << " to persistent.");
                ao.get()->pushToPersistent(_subEventWeights);
            }
            MSG_TRACE("AnalysisHandler::analyze(): finished pushing " << a->name() << "'s objects to persistent.");
        }

        _eventNumber = ge.event_number();

        MSG_DEBUG("nominal event # " << _eventCounter.get()->_persistent[0]->numEntries());
        MSG_DEBUG("nominal sum of weights: " << _eventCounter.get()->_persistent[0]->sumW());
        MSG_DEBUG("Event has " << _subEventWeights.size() << " sub events.");
        _subEventWeights.clear();
    }


    MSG_TRACE("starting new sub event");
    _eventCounter.get()->newSubEvent();

    for (const AnaHandle& a : _analyses) {
        for (auto ao : a->analysisObjects()) {
            ao.get()->newSubEvent();
        }
    }

    _subEventWeights.push_back(event.weights());
    MSG_DEBUG("Analyzing subevent #" << _subEventWeights.size() - 1 << ".");

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

    if ( _dumpPeriod > 0 && numEvents()%_dumpPeriod == 0 ) {
      MSG_INFO("Dumping intermediate results to " << _dumpFile << ".");
      _dumping = true;
      finalize();
      _dumping = false;
      writeData(_dumpFile);
    }

  }


  void AnalysisHandler::analyze(const GenEvent* ge) {
    if (ge == nullptr) {
      MSG_ERROR("AnalysisHandler received null pointer to GenEvent");
      //throw Error("AnalysisHandler received null pointer to GenEvent");
    }
    analyze(*ge);
  }


  void AnalysisHandler::finalize() {
    if (!_initialised) return;
    MSG_INFO("Finalising analyses");

    // First push all analyses' objects to persistent
    MSG_TRACE("AnalysisHandler::finalize(): Pushing analysis objects to persistent.");
    _eventCounter.get()->pushToPersistent(_subEventWeights);
    for (const AnaHandle& a : _analyses) {
      for (auto ao : a->analysisObjects())
        ao.get()->pushToPersistent(_subEventWeights);
    }

    // First we make copies of all analysis objects.
    map<string,AnalysisObjectPtr> backupAOs;
    for (auto ao : getYodaAOs(false, true, false) )
      backupAOs[ao->path()] = AnalysisObjectPtr(ao->newclone());

    // Finalize all the histograms
    for (const AnaHandle& a : _analyses) {
      // a->setCrossSection(_xs);
      for (size_t iW = 0; iW < numWeights(); iW++) {
        _eventCounter.get()->setActiveWeightIdx(iW);
        _xs.get()->setActiveWeightIdx(iW);
        for (auto ao : a->analysisObjects())
          ao.get()->setActiveWeightIdx(iW);

        MSG_TRACE("Running " << a->name() << "::finalize() for weight " << iW << ".");

        try {
          if ( !_dumping || a->info().reentrant() )  a->finalize();
          else if ( _dumping == 1 && iW == 0 )
            MSG_INFO("Skipping periodic dump of " << a->name()
                     << " as it is not declared reentrant.");
        } catch (const Error& err) {
          cerr << "Error in " << a->name() << "::finalize method: " << err.what() << endl;
          exit(1);
        }

      }
    }

    // Now we copy all analysis objects to the list of finalized
    // ones, and restore the value to their original ones.
    _finalizedAOs.clear();
    for ( auto ao : getYodaAOs(false, false, false) )
      _finalizedAOs.push_back(AnalysisObjectPtr(ao->newclone()));
    for ( auto ao : getYodaAOs(false, true, false) ) {
      // TODO: This should be possible to do in a nicer way, with a flag etc.
      if (ao->path().find("/FINAL") != std::string::npos) continue;
      auto aoit = backupAOs.find(ao->path());
      if ( aoit == backupAOs.end() ) {
        AnaHandle ana = analysis(split(ao->path(), "/")[0]);
        if ( ana ) ana->removeAnalysisObject(ao->path());
      } else
        copyao(aoit->second, ao);
    }

    // Print out number of events processed
    const int nevts = numEvents();
    MSG_INFO("Processed " << nevts << " event" << (nevts != 1 ? "s" : ""));

    // // Delete analyses
    // MSG_DEBUG("Deleting analyses");
    // _analyses.clear();

    // Print out MCnet boilerplate
    cout << '\n';
    cout << "The MCnet usage guidelines apply to Rivet: see http://www.montecarlonet.org/GUIDELINES" << '\n';
    cout << "Please acknowledge plots made with Rivet analyses, and cite arXiv:1003.0694 (http://arxiv.org/abs/1003.0694)" << endl;
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
          MSG_WARNING("Error in option specification. Skipping analysis "
                      << analysisname);
          return *this;
        }
        if ( !analysis->info().validOption(opt[0], opt[1]) ) {
          MSG_WARNING("Cannot set option '" << opt[0] << "' to '" << opt[1]
                      << "'. Skipping analysis " << analysisname);
          return *this;
        }
        opts[opt[0]] = opt[1];
      }
      for ( auto opt: opts) {
        analysis->_options[opt.first] = opt.second;
        analysis->_optstring += ":" + opt.first + "=" + opt.second;
      }
      for (const AnaHandle& a : _analyses) {
        if (a->name() == analysis->name() ) {
          MSG_WARNING("Analysis '" << analysisname << "' already registered: skipping duplicate");
          return *this;
        }
      }
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


  void AnalysisHandler::addData(const std::vector<YODA::AnalysisObjectPtr>& aos) {
    for (const YODA::AnalysisObjectPtr ao : aos) {
      string path = ao->path();
      if ( path.substr(0, 5) != "/RAW/" ) {
        _orphanedPreloads.push_back(ao);
        continue;
      }

      path = path.substr(4);
      ao->setPath(path);
      if (path.size() > 1) { // path > "/"
        try {
          const string ananame =  ::split(path, "/")[0];
          AnaHandle a = analysis(ananame);
          //MultiweightAOPtr mao = ????; /// @todo generate right Multiweight object from ao
          a->addAnalysisObject(mao); /// @todo Need to statistically merge...
        } catch (const Error& e) {
          MSG_TRACE("Adding analysis object " << path <<
                    " to the list of orphans.");
          _orphanedPreloads.push_back(ao);
        }
      }
    }
  }

  void AnalysisHandler::stripOptions(AnalysisObjectPtr ao,
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




  void AnalysisHandler::
  mergeYodas(const vector<string> & aofiles, const vector<string> & delopts, bool equiv) {
    vector< vector<AnalysisObjectPtr> > aosv;
    vector<double> xsecs;
    vector<double> xsecerrs;
    vector<CounterPtr> sows;
    set<string> ananames;
     _eventcounter.reset();

    // First scan all files and extract analysis objects and add the
    // corresponding anayses..
    for ( auto file : aofiles ) {
      Scatter1DPtr xsec;
      CounterPtr sow;

      // For each file make sure that cross section and sum-of-weights
      // objects are present and stor all RAW ones in a vector;
      vector<AnalysisObjectPtr> aos;
      try {
        /// @todo Use new YODA SFINAE to fill the smart ptr vector directly
        vector<YODA::AnalysisObject*> aos_raw;
        YODA::read(file, aos_raw);
        for (AnalysisObject* aor : aos_raw) {
          AnalysisObjectPtr ao = AnalysisObjectPtr(aor);
          if ( ao->path().substr(0, 5) != "/RAW/" ) continue;
          ao->setPath(ao->path().substr(4));
          if ( ao->path() == "/_XSEC" )
            xsec = dynamic_pointer_cast<Scatter1D>(ao);
          else if ( ao->path() == "/_EVTCOUNT" )
            sow = dynamic_pointer_cast<Counter>(ao);
          else {
            stripOptions(ao, delopts);
            string ananame = split(ao->path(), "/")[0];
            if ( ananames.insert(ananame).second ) addAnalysis(ananame);
            aos.push_back(ao);
          }
        }
        if ( !xsec || !sow ) {
          MSG_ERROR( "Error in AnalysisHandler::mergeYodas: The file " << file
                     << " did not contain weights and cross section info.");
          exit(1);
        }
        xsecs.push_back(xsec->point(0).x());
        sows.push_back(sow);
	xsecerrs.push_back(sqr(xsec->point(0).xErrAvg()));
        _eventcounter += *sow;
        sows.push_back(sow);
        aosv.push_back(aos);
      } catch (...) { //< YODA::ReadError&
        throw UserError("Unexpected error in reading file: " + file);
      }
    }

    // Now calculate the scale to be applied for all bins in a file
    // and get the common cross section and sum of weights.
    _xs = _xserr = 0.0;
    for ( int i = 0, N = sows.size(); i < N; ++i ) {
      double effnent = sows[i]->effNumEntries();
      _xs += (equiv? effnent: 1.0)*xsecs[i];
      _xserr += (equiv? sqr(effnent): 1.0)*xsecerrs[i];
    }

    vector<double> scales(sows.size(), 1.0);
    if ( equiv ) {
      _xs /= _eventcounter.effNumEntries();
      _xserr = sqrt(_xserr)/_eventcounter.effNumEntries();
    } else {
      _xserr = sqrt(_xserr);
      for ( int i = 0, N = sows.size(); i < N; ++i )
        scales[i] = (_eventcounter.sumW()/sows[i]->sumW())*(xsecs[i]/_xs);
    }

    // Initialize the analyses allowing them to book analysis objects.
    for (AnaHandle a : _analyses) {
      MSG_DEBUG("Initialising analysis: " << a->name());
      if ( !a->info().reentrant() )
        MSG_WARNING("Analysis " << a->name() << " has not been validated to have "
                    << "a reentrant finalize method. The result is unpredictable.");
      try {
        // Allow projection registration in the init phase onwards
        a->_allowProjReg = true;
        cerr << "sqrtS " << sqrtS() << endl;
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
    // Get a list of all anaysis objects to handle.
    map<string,AnalysisObjectPtr> current;
    for ( auto ao : getData(false, true) ) current[ao->path()] = ao;
    // Go through all objects to be merged and add them to current
    // after appropriate scaling.
    for ( int i = 0, N = aosv.size(); i < N; ++i)
      for ( auto ao : aosv[i] ) {
        if ( ao->path() == "/_XSEC" || ao->path() == "_EVTCOUNT" ) continue;
	auto aoit = current.find(ao->path());
        if ( aoit == current.end() ) {
          MSG_WARNING("" << ao->path() << " was not properly booked.");
          continue;
        }
        if ( !addaos(aoit->second, ao, scales[i]) )
          MSG_WARNING("Cannot merge objects with path " << ao->path()
                      <<" of type " << ao->annotation("Type") );
      }
    // Now we can simply finalize() the analysis, leaving the
    // controlling program to write it out some yoda-file.
    finalize();

  }


  void AnalysisHandler::readData(const string& filename) {
    vector<YODA::AnalysisObjectPtr> aos;
    try {
      /// @todo Use new YODA SFINAE to fill the smart ptr vector directly
      vector<YODA::AnalysisObject*> aos_raw;
      YODA::ReaderYODA::read(filename, aos_raw);
      for (YODA::AnalysisObject* aor : aos_raw) aos.push_back(YODA::AnalysisObjectPtr(aor));
    //} catch (const YODA::ReadError & e) {
    } catch (...) { //< YODA::ReadError&
      throw UserError("Unexpected error in reading file: " + filename);
    }
    if (!aos.empty()) addData(aos);
  }


  vector<MultiweightAOPtr> AnalysisHandler::getRivetAOs() const {
      vector<MultiweightAOPtr> rtn;

      for (AnaHandle a : _analyses) {
          for (const auto & ao : a->analysisObjects()) {
              rtn.push_back(ao);
          }
      }

      rtn.push_back(_eventCounter);
      rtn.push_back(_xs);

      return rtn;
  }

  vector<YODA::AnalysisObjectPtr> AnalysisHandler::getYodaAOs(bool includeorphans,
                                                              bool includetmps,
                                                              bool usefinalized) const {
      vector<YODA::AnalysisObjectPtr> rtn;
      if (usefinalized)
        rtn = _finalizedAOs;
      else {
        for (auto rao : getRivetAOs()) {
          // need to set the index
          // before we can search the PATH
          rao.get()->setActiveWeightIdx(_defaultWeightIdx);
          // Exclude paths from final write-out if they contain a "TMP" layer (i.e. matching "/TMP/")
          if (!includetmps && rao->path().find("/TMP/") != string::npos)
            continue;

          for (size_t iW = 0; iW < numWeights(); iW++) {
            rao.get()->setActiveWeightIdx(iW);
            rtn.push_back(rao.get()->activeYODAPtr());
          }
        }
      }

      // Sort histograms alphanumerically by path before write-out
      sort(rtn.begin(), rtn.end(),
           [](YODA::AnalysisObjectPtr a, YODA::AnalysisObjectPtr b) {
                return a->path() < b->path();
            }
          );

      return rtn;
  }


  vector<YODA::AnalysisObjectPtr> AnalysisHandler::getData(bool includeorphans,
                                                           bool includetmps,
                                                           bool usefinalized) const {
    return getYodaAOs(includeorphans, includetmps, usefinalized);
  }


  void AnalysisHandler::writeData(const string& filename) const {
    vector<YODA::AnalysisObjectPtr> out = _finalizedAOs;
    set<string> finalana;
    for ( auto ao : out) finalana.insert(ao->path());
    out.reserve(2*out.size());
    vector<YODA::AnalysisObjectPtr> aos = getData(false, true, false);

    if ( _dumping ) {
      for ( auto ao : aos ) {
        if ( finalana.find(ao->path()) == finalana.end() )
          out.push_back(AnalysisObjectPtr(ao->newclone()));
      }
    }

    for ( auto ao : aos ) {
      ao = YODA::AnalysisObjectPtr(ao->newclone());
      ao->setPath("/RAW" + ao->path());
      out.push_back(ao);
    }

    try {
      YODA::WriterYODA::write(filename, aos.begin(), aos.end());
    } catch (...) { //< YODA::WriteError&
      throw UserError("Unexpected error in writing file: " + filename);
    }
  }


  string AnalysisHandler::runName() const { return _runname; }
  size_t AnalysisHandler::numEvents() const { return _eventCounter->numEntries(); }


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


  AnalysisHandler& AnalysisHandler::setCrossSection(double xs, double xserr) {
    _xs = Scatter1DPtr(weightNames(), Scatter1D("_XSEC"));
    _eventCounter.get()->setActiveWeightIdx(_defaultWeightIdx);
    double nomwgt = sumOfWeights();

    // The cross section of each weight variation is the nominal cross section
    // times the sumOfWeights(variation) / sumOfWeights(nominal).
    // This way the cross section will work correctly
    for (size_t iW = 0; iW < numWeights(); iW++) {
      _eventCounter.get()->setActiveWeightIdx(iW);
      double s = sumOfWeights() / nomwgt;
      _xs.get()->setActiveWeightIdx(iW);
      _xs->addPoint(xs*s, xserr*s);
    }

    _eventCounter.get()->unsetActiveWeight();
    _xs.get()->unsetActiveWeight();
    return *this;
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
