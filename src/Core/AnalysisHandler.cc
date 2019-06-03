// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/ParticleName.hh"
#include "Rivet/Tools/BeamConstraint.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Beam.hh"
#include "YODA/IO.h"
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
      _eventcounter("/_EVTCOUNT"),
      _xs(NAN), _xserr(NAN),
      _initialised(false), _ignoreBeams(false), _dumpPeriod(0), _dumping(false),
      _defaultWeightIdx(0)
  {  }


  AnalysisHandler::~AnalysisHandler()
  {  }


  Log& AnalysisHandler::getLog() const {
    return Log::getLog("Rivet.Analysis.Handler");
  }


  /// http://stackoverflow.com/questions/4654636/how-to-determine-if-a-string-is-a-number-with-c
  namespace {
    bool is_number(const std::string& s) {
      std::string::const_iterator it = s.begin();
      while (it != s.end() && std::isdigit(*it)) ++it;
      return !s.empty() && it == s.end();
    }
  }

  /// Check if any of the weightnames is not a number
  bool AnalysisHandler::haveNamedWeights() const {
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
    for (const AnaHandle a : analyses()) {
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
    bool strip = ( getEnvParam("RIVET_STRIP_HEPMC", string("NOOOO") ) != "NOOOO" );
    Event event(ge, strip);

    // Weights
    /// @todo Drop this / just report first weight when we support multiweight events
    _eventcounter.fill(event.weight());
    MSG_DEBUG("Event #" << _eventcounter.numEntries() << " weight = " << event.weight());

    // Cross-section
    #if defined ENABLE_HEPMC_3
    if (ge.cross_section()) {
      //@todo HepMC3::GenCrossSection methods aren't const accessible :(
      RivetHepMC::GenCrossSection gcs = *(event.genEvent()->cross_section());
      _xs = gcs.xsec();
      _xserr = gcs.xsec_err();
    }
    #elif defined HEPMC_HAS_CROSS_SECTION
    if (ge.cross_section()) {
      _xs = ge.cross_section()->cross_section();
      _xserr = ge.cross_section()->cross_section_error();
    }
    #endif

    // Won't happen for first event because _eventNumber is set in init()
    if (_eventNumber != ge.event_number()) {
        /// @todo Can we get away with not passing a matrix?

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

    if ( _dumpPeriod > 0 && numEvents()%_dumpPeriod == 0 ) {
      MSG_INFO("Dumping intermediate results to " << _dumpFile << ".");
      _dumping = numEvents()/_dumpPeriod;
      finalize();
      writeData(_dumpFile);
      _dumping = 0;
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
    for (auto ao : getData(false, true, false) )
      backupAOs[ao->path()] = AnalysisObjectPtr(ao->newclone());

    // Now we run the (re-entrant) finalize() functions for all analyses.
    MSG_INFO("Finalising analyses");
    for (AnaHandle a : analyses()) {
      a->setCrossSection(_xs);
      try {
        if ( !_dumping || a->info().reentrant() )  a->finalize();
        else if ( _dumping == 1 )
          MSG_INFO("Skipping finalize in periodic dump of " << a->name()
                   << " as it is not declared reentrant.");
      } catch (const Error& err) {
        cerr << "Error in " << a->name() << "::finalize method: " << err.what() << endl;
        exit(1);
      }
    }

    // Now we copy all analysis objects to the list of finalized
    // ones, and restore the value to their original ones.
    _finalizedAOs.clear();
    for ( auto ao : getYodaAOs(false, false, false) )
      _finalizedAOs.push_back(YODA::AnalysisObjectPtr(ao->newclone()));
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
    cout << endl;
    cout << "The MCnet usage guidelines apply to Rivet: see http://www.montecarlonet.org/GUIDELINES" << endl;
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
          MSG_WARNING("Error in option specification. Skipping analysis " << analysisname);
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
          /// @todo FIXXXXX
          //MultiweightAOPtr mao = ????; /// @todo generate right Multiweight object from ao
          //a->addAnalysisObject(mao); /// @todo Need to statistically merge...
        } catch (const Error& e) {
          MSG_TRACE("Adding analysis object " << path <<
                    " to the list of orphans.");
          _orphanedPreloads.push_back(ao);
        }
      }
    }
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


  void AnalysisHandler::mergeYodas(const vector<string> & ,
                                   const vector<string> & , bool ) {}
  // void AnalysisHandler::mergeYodas(const vector<string> & aofiles,
  //                                  const vector<string> & delopts, bool equiv) {
  //   vector< vector<YODA::AnalysisObjectPtr> > aosv;
  //   vector<double> xsecs;
  //   vector<double> xsecerrs;
  //   vector<CounterPtr> sows;
  //   set<string> ananames;
  //   _eventCounter->reset();

  //   // First scan all files and extract analysis objects and add the
  //   // corresponding analyses..
  //   for (const string& file : aofiles ) {
  //     Scatter1DPtr xsec;
  //     CounterPtr sow;

  //     // For each file make sure that cross section and sum-of-weights
  //     // objects are present and stor all RAW ones in a vector;
  //     vector<YODA::AnalysisObjectPtr> aos;
  //     try {
  //       /// @todo Use new YODA SFINAE to fill the smart ptr vector directly
  //       vector<YODA::AnalysisObject*> aos_raw;
  //       YODA::read(file, aos_raw);
  //       for (YODA::AnalysisObject* aor : aos_raw) {
  //         YODA::AnalysisObjectPtr ao(aor);
  //         if ( ao->path().substr(0, 5) != "/RAW/" ) continue;
  //         ao->setPath(ao->path().substr(4));
  //         if ( ao->path() == "/_XSEC" )
  //           xsec = dynamic_pointer_cast<Scatter1D>(ao);
  //         else if ( ao->path() == "/_EVTCOUNT" )
  //           sow = dynamic_pointer_cast<Counter>(ao);
  //         else {
  //           stripOptions(ao, delopts);
  //           string ananame = split(ao->path(), "/")[0];
  //           if ( ananames.insert(ananame).second ) addAnalysis(ananame);
  //           aos.push_back(ao);
  //         }
  //       }
  //       if ( !xsec || !sow ) {
  //         MSG_ERROR( "Error in AnalysisHandler::mergeYodas: The file " << file
  //                    << " did not contain weights and cross section info.");
  //         exit(1);
  //       }
  //       xsecs.push_back(xsec->point(0).x());
  //       sows.push_back(sow);
  //       xsecerrs.push_back(sqr(xsec->point(0).xErrAvg()));
  //       _eventCounter->operator+=(*sow); //< HAHAHAHAHA!
  //       sows.push_back(sow);
  //       aosv.push_back(aos);
  //     } catch (...) { //< YODA::ReadError&
  //       throw UserError("Unexpected error in reading file: " + file);
  //     }
  //   }

  //   // Now calculate the scale to be applied for all bins in a file
  //   // and get the common cross section and sum of weights.
  //   _xs = _xserr = 0.0;
  //   for ( int i = 0, N = sows.size(); i < N; ++i ) {
  //     double effnent = sows[i]->effNumEntries();
  //     _xs += (equiv? effnent: 1.0)*xsecs[i];
  //     _xserr += (equiv? sqr(effnent): 1.0)*xsecerrs[i];
  //   }

  //   vector<double> scales(sows.size(), 1.0);
  //   if ( equiv ) {
  //     _xs /= _eventCounter.effNumEntries();
  //     _xserr = sqrt(_xserr)/_eventCounter.effNumEntries();
  //   } else {
  //     _xserr = sqrt(_xserr);
  //     for ( int i = 0, N = sows.size(); i < N; ++i )
  //       scales[i] = (_eventCounter.sumW()/sows[i]->sumW())*(xsecs[i]/_xs);
  //   }

  //   // Initialize the analyses allowing them to book analysis objects.
  //   for (AnaHandle a : _analyses) {
  //     MSG_DEBUG("Initialising analysis: " << a->name());
  //     if ( !a->info().reentrant() )
  //       MSG_WARNING("Analysis " << a->name() << " has not been validated to have "
  //                   << "a reentrant finalize method. The result is unpredictable.");
  //     try {
  //       // Allow projection registration in the init phase onwards
  //       a->_allowProjReg = true;
  //       cerr << "sqrtS " << sqrtS() << endl;
  //       a->init();
  //       //MSG_DEBUG("Checking consistency of analysis: " << a->name());
  //       //a->checkConsistency();
  //     } catch (const Error& err) {
  //       cerr << "Error in " << a->name() << "::init method: " << err.what() << endl;
  //       exit(1);
  //     }
  //     MSG_DEBUG("Done initialising analysis: " << a->name());
  //   }
  //   _initialised = true;
  //   // Get a list of all anaysis objects to handle.
  //   map<string,AnalysisObjectPtr> current;
  //   for ( auto ao : getData(false, true) ) current[ao->path()] = ao;
  //   // Go through all objects to be merged and add them to current
  //   // after appropriate scaling.
  //   for ( int i = 0, N = aosv.size(); i < N; ++i)
  //     for ( auto ao : aosv[i] ) {
  //       if ( ao->path() == "/_XSEC" || ao->path() == "_EVTCOUNT" ) continue;
  //       auto aoit = current.find(ao->path());
  //       if ( aoit == current.end() ) {
  //         MSG_WARNING("" << ao->path() << " was not properly booked.");
  //         continue;
  //       }
  //       if ( !addaos(aoit->second, ao, scales[i]) )
  //         MSG_WARNING("Cannot merge objects with path " << ao->path()
  //                     <<" of type " << ao->annotation("Type") );
  //     }
  //   // Now we can simply finalize() the analysis, leaving the
  //   // controlling program to write it out some yoda-file.
  //   finalize();
  // }


  void AnalysisHandler::readData(const string& filename) {
    vector<YODA::AnalysisObjectPtr> aos;
    try {
      /// @todo Use new YODA SFINAE to fill the smart ptr vector directly
      vector<YODA::AnalysisObject*> aos_raw;
      YODA::read(filename, aos_raw);
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
    vector<AnalysisObjectPtr> aos = getData(false, true);
    for ( auto ao : aos ) {
      ao = YODA::AnalysisObjectPtr(ao->newclone());
      ao->setPath("/RAW" + ao->path());
      out.push_back(ao);
    }

    try {
      YODA::write(filename, aos.begin(), aos.end());
    } catch (...) { //< YODA::WriteError&
      throw UserError("Unexpected error in writing file: " + filename);
    }
  }


  string AnalysisHandler::runName() const { return _runname; }
  size_t AnalysisHandler::numEvents() const { return _eventCounter->numEntries(); }


  std::vector<std::string> AnalysisHandler::analysisNames() const {
    std::vector<std::string> rtn;
    for (AnaHandle a : analyses()) {
      rtn.push_back(a->name());
    }
    return rtn;
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
    double nomwgt = sumW();

    // The cross section of each weight variation is the nominal cross section
    // times the sumW(variation) / sumW(nominal).
    // This way the cross section will work correctly
    for (size_t iW = 0; iW < numWeights(); iW++) {
      _eventCounter.get()->setActiveWeightIdx(iW);
      double s = sumW() / nomwgt;
      _xs.get()->setActiveWeightIdx(iW);
      _xs->addPoint(xs*s, xserr*s);
    }

    _eventCounter.get()->unsetActiveWeight();
    _xs.get()->unsetActiveWeight();
    return *this;
  }


  bool AnalysisHandler::hasCrossSection() const {
    return (!std::isnan(crossSection()));
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


}
