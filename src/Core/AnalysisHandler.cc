// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/ParticleName.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projections/Beam.hh"
#include "LWH/AIManagedObject.h"

using namespace AIDA;

namespace Rivet {


  AnalysisHandler::AnalysisHandler(const string& runname)
    : _runname(runname), _numEvents(0),
      _sumOfWeights(0.0), _xs(-1.0),
      _initialised(false)
  {
    _theAnalysisFactory.reset( createAnalysisFactory() );
    _setupFactories();
  }


  AnalysisHandler::AnalysisHandler(const string& basefilename,
                                   const string& runname, HistoFormat storetype)
    : _runname(runname), _numEvents(0),
      _sumOfWeights(0.0), _xs(-1.0),
      _initialised(false)
  {
    cerr << "AnalysisHandler(basefilename, runname, format) constructor is deprecated: "
         << "please migrate your code to use the one-arg constructor" << endl;
    _theAnalysisFactory.reset( createAnalysisFactory() );
    _setupFactories(basefilename, storetype);
  }


  // AnalysisHandler::AnalysisHandler(IAnalysisFactory& afac, string basefilename,
  //                                  string runname, HistoFormat storetype)
  //   : _runname(runname), _numEvents(0),
  //     _sumOfWeights(0.0), _xs(-1.0), _initialised(false),
  //     _theAnalysisFactory(&afac)
  // {
  //   _setupFactories(basefilename, storetype);
  //   initializeParticleNames();
  // }


  AnalysisHandler::~AnalysisHandler()
  {  }


  Log& AnalysisHandler::getLog() {
    return Log::getLog("Rivet.Analysis.Handler");
  }


  void AnalysisHandler::init(const GenEvent& ge) {
    assert(!_initialised);
    setRunBeams(Rivet::beams(ge));
    getLog() << Log::DEBUG << "Initialising the analysis handler" << endl;
    _numEvents = 0;
    _sumOfWeights = 0.0;

    // Check that analyses are beam-compatible
    const size_t num_anas_requested = analysisNames().size();
    removeIncompatibleAnalyses(beamIds());
    foreach (const AnaHandle a, analyses()) {
      if (toUpper(a->status()) != "VALIDATED") {
        getLog() << Log::WARN
                 << "Analysis '" << a->name() << "' is unvalidated: be careful!" << endl;
      }
    }
    if (num_anas_requested > 0 && analysisNames().size() == 0) {
      getLog() << Log::ERROR
               << "All analyses were incompatible with the first event's beams\n"
               << "Exiting, since this probably isn't intentional!" << endl;
      exit(1);
    }

    foreach (AnaHandle a, _analyses) {
      getLog() << Log::DEBUG << "Initialising analysis: " << a->name() << endl;
      try {
        // Allow projection registration in the init phase onwards
        a->_allowProjReg = true;
        a->init();
        //getLog() << Log::DEBUG << "Checking consistency of analysis: " << a->name() << endl;
        //a->checkConsistency();
      } catch (const Error& err) {
        getLog() << Log::ERROR << "Error in " << a->name() << "::init method: "
                 << err.what() << endl;
        exit(1);
      }
      getLog() << Log::DEBUG << "Done initialising analysis: " << a->name() << endl;
    }
    _initialised = true;
    getLog() << Log::DEBUG << "Analysis handler initialised" << endl;
  }


  void AnalysisHandler::analyze(const GenEvent& ge) {
    // Call init with event as template if not already initialised
    if (!_initialised) {
      init(ge);
    }
    // Proceed with event analysis
    assert(_initialised);
    // Ensure that beam details match those from first event
    const PdgIdPair beams = Rivet::beamIds(ge);
    const double sqrts = Rivet::sqrtS(ge);
    if (!compatible(beams, _beams) || !fuzzyEquals(sqrts, sqrtS())) {
      getLog() << Log::ERROR << "Event beams mismatch: "
               << toBeamsString(beams) << " @ " << sqrts/GeV << " GeV" << " vs. first beams "
               << this->beams() << " @ " << this->sqrtS()/GeV << " GeV" << endl;
      exit(1);
    }


    Event event(ge);
    _numEvents++;
    // Weights
    const double weight = event.weight();
    _sumOfWeights += weight;
    getLog() << Log::DEBUG << "Event #" << _numEvents << " weight = " << weight << endl;
    #ifdef HEPMC_HAS_CROSS_SECTION
    if (ge.cross_section()) {
      const double xs = ge.cross_section()->cross_section();
      setCrossSection(xs);
    }
    #endif
    foreach (AnaHandle a, _analyses) {
      //getLog() << Log::DEBUG << "About to run analysis " << a->name() << endl;
      try {
        a->analyze(event);
      } catch (const Error& err) {
        getLog() << Log::ERROR << "Error in " << a->name() << "::analyze method: "
                 << err.what() << endl;
        exit(1);
      }
      //getLog() << Log::DEBUG << "Finished running analysis " << a->name() << endl;
    }
  }


  void AnalysisHandler::finalize() {
    assert(_initialised);
    getLog() << Log::INFO << "Finalising analyses" << endl;
    foreach (AnaHandle a, _analyses) {
      try {
        a->finalize();
      } catch (const Error& err) {
        getLog() << Log::ERROR << "Error in " << a->name() << "::finalize method: "
                 << err.what() << endl;
        exit(1);
      }
    }

    // Print out number of events processed
    getLog() << Log::INFO << "Processed " << _numEvents << " event" << (_numEvents == 1 ? "" : "s") << endl;

    // Change AIDA histos into data point sets
    getLog() << Log::DEBUG << "Converting histograms to scatter plots" << endl;
    assert(_theTree != 0);
    _normalizeTree(tree());

    // Delete analyses
    getLog() << Log::DEBUG << "Deleting analyses" << endl;
    _analyses.clear();

    // Print out MCnet boilerplate
    cout << endl;
    cout << "The MCnet usage guidelines apply to Rivet: see http://www.montecarlonet.org/GUIDELINES" << endl;
    cout << "Please acknowledge plots made with Rivet analyses, and cite arXiv:1003.0694 (http://arxiv.org/abs/1003.0694)" << endl;
  }


  AnalysisHandler& AnalysisHandler::addAnalysis(const string& analysisname) {
    // Check for a duplicate analysis
    /// @todo Might we want to be able to run an analysis twice, with different params?
    ///       Requires avoiding histo tree clashes, i.e. storing the histos on the analysis objects.
    foreach (const AnaHandle& a, _analyses) {
      if (a->name() == analysisname) {
        getLog() << Log::WARNING 
		 << "Analysis '" << analysisname 
		 << "' already registered: skipping duplicate" << endl;
        return *this;
      }
    }
    AnaHandle analysis( AnalysisLoader::getAnalysis(analysisname) );
    if (analysis.get() != 0) { // < Check for null analysis.
      getLog() << Log::DEBUG << "Adding analysis '" << analysisname << "'" << endl;
      analysis->_analysishandler = this;
      _analyses.insert(analysis);
    }
    else {
      getLog() << Log::WARNING 
	       << "Analysis '" << analysisname 
	       << "' not found." << endl;
    }
    return *this;
  }


  AnalysisHandler& AnalysisHandler::removeAnalysis(const string& analysisname) {
    shared_ptr<Analysis> toremove;
    foreach (const AnaHandle a, _analyses) {
      if (a->name() == analysisname) {
        toremove.reset( a.get() );
        break;
      }
    }
    if (toremove.get() != 0) {
      getLog() << Log::DEBUG << "Removing analysis '" << analysisname << "'" << endl;
      _analyses.erase(toremove);
    }
    return *this;
  }


  AnalysisHandler& AnalysisHandler::removeIncompatibleAnalyses(const PdgIdPair& beams) {
    vector<string> anamestodelete;
    foreach (const AnaHandle a, _analyses) {
      if (! a->isCompatible(beams)) {
        anamestodelete.push_back(a->name());
      }
    }
    foreach (const string& aname, anamestodelete) {
      getLog() << Log::WARN << "Removing incompatible analysis '"
               << aname << "'" << endl;
      removeAnalysis(aname);
    }
    return *this;
  }


  void AnalysisHandler::_setupFactories(const string& basefilename, HistoFormat storetype) {
    string filename(basefilename), storetypestr("");
    if (storetype == AIDAML) {
      if (!endsWith(filename, ".aida")) filename += ".aida";
      storetypestr = "xml";
    } else if (storetype == FLAT) {
      if (!endsWith(filename, ".data")) filename += ".data";
      storetypestr = "flat";
    } else if (storetype == ROOT) {
      if (!endsWith(filename, ".root")) filename += ".root";
      storetypestr = "root";
    }
    _theTreeFactory = _theAnalysisFactory->createTreeFactory();
    _theTree = _theTreeFactory->create(filename, storetypestr, false, true);
    _theHistogramFactory = _theAnalysisFactory->createHistogramFactory(tree());
    _theDataPointSetFactory = _theAnalysisFactory->createDataPointSetFactory(tree());
  }


  void AnalysisHandler::_setupFactories() {
    _theTreeFactory = _theAnalysisFactory->createTreeFactory();
    _theTree = _theTreeFactory->create();
    _theHistogramFactory = _theAnalysisFactory->createHistogramFactory(tree());
    _theDataPointSetFactory = _theAnalysisFactory->createDataPointSetFactory(tree());
  }


  void AnalysisHandler::commitData() {
    tree().commit();
  }


  void AnalysisHandler::writeData(const string& filename) {
    tree().commit(filename);
  }


  void AnalysisHandler::_normalizeTree(ITree& tree) {
    Log& log = getLog();
    const vector<string> paths = tree.listObjectNames("/", true); // args set recursive listing
    log << Log::TRACE << "Number of objects in AIDA tree = " << paths.size() << endl;
    const string tmpdir = "/RivetNormalizeTmp";
    tree.mkdir(tmpdir);
    foreach (const string& path, paths) {

      IManagedObject* hobj = tree.find(path);
      if (hobj) {

        // Weird seg fault on SLC4 when trying to dyn cast an IProfile ptr to a IHistogram
        // Fix by attempting to cast to IProfile first, only try IHistogram if it fails.
        IHistogram1D* histo = 0;
        IProfile1D* prof = dynamic_cast<IProfile1D*>(hobj);
        if (!prof) histo = dynamic_cast<IHistogram1D*>(hobj);

        // If it's a normal histo:
        if (histo) {
          log << Log::TRACE << "Converting histo " << path << " to DPS" << endl;
          tree.mv(path, tmpdir);
          const size_t lastslash = path.find_last_of("/");
          const string basename = path.substr(lastslash+1, path.length() - (lastslash+1));
          const string tmppath = tmpdir + "/" + basename;
          IHistogram1D* tmphisto = dynamic_cast<IHistogram1D*>(tree.find(tmppath));
          if (tmphisto) {
            //getLog() << Log::TRACE << "Temp histo " << tmppath << " exists" << endl;
            datapointsetFactory().create(path, *tmphisto);
          }
          tree.rm(tmppath);
        }
        // If it's a profile histo:
        else if (prof) {
          log << Log::TRACE << "Converting profile histo " << path << " to DPS" << endl;
          tree.mv(path, tmpdir);
          const size_t lastslash = path.find_last_of("/");
          const string basename = path.substr(lastslash+1, path.length() - (lastslash+1));
          const string tmppath = tmpdir + "/" + basename;
          IProfile1D* tmpprof = dynamic_cast<IProfile1D*>(tree.find(tmppath));
          if (tmpprof) {
            //getLog() << Log::TRACE << "Temp profile histo " << tmppath << " exists" << endl;
            datapointsetFactory().create(path, *tmpprof);
          }
          tree.rm(tmppath);
        }

      }

    }
    tree.rmdir(tmpdir);
  }


  string AnalysisHandler::runName() const { return _runname; }
  size_t AnalysisHandler::numEvents() const { return _numEvents; }
  double AnalysisHandler::sumOfWeights() const { return _sumOfWeights; }


  void AnalysisHandler::setSumOfWeights(const double& sum) {
    _sumOfWeights=sum;
  }


  std::vector<std::string> AnalysisHandler::analysisNames() const {
    std::vector<std::string> rtn;
    foreach (AnaHandle a, _analyses) {
      rtn.push_back(a->name());
    }
    return rtn;
  }


  AnalysisHandler& AnalysisHandler::addAnalyses(const std::vector<std::string>& analysisnames) {
    foreach (const string& aname, analysisnames) {
      //getLog() << Log::DEBUG << "Adding analysis '" << aname << "'" << endl;
      addAnalysis(aname);
    }
    return *this;
  }


  AnalysisHandler& AnalysisHandler::removeAnalyses(const std::vector<std::string>& analysisnames) {
    foreach (const string& aname, analysisnames) {
      removeAnalysis(aname);
    }
    return *this;
  }



  AIDA::IAnalysisFactory& AnalysisHandler::analysisFactory() {
    return *_theAnalysisFactory;
  }


  AIDA::ITree& AnalysisHandler::tree() {
    return *_theTree;
  }


  AIDA::IHistogramFactory& AnalysisHandler::histogramFactory() {
    return *_theHistogramFactory;
  }


  AIDA::IDataPointSetFactory& AnalysisHandler::datapointsetFactory() {
    return *_theDataPointSetFactory;
  }


  bool AnalysisHandler::needCrossSection() const {
    bool rtn = false;
    foreach (const AnaHandle a, _analyses) {
      if (!rtn) rtn = a->needsCrossSection();
      if (rtn) break;
    }
    return rtn;
  }


  AnalysisHandler& AnalysisHandler::setCrossSection(double xs) {
    _xs = xs;
    foreach (AnaHandle a, _analyses) {
      a->setCrossSection(xs);
    }
    return *this;
  }


  bool AnalysisHandler::hasCrossSection() const {
    return (crossSection() >= 0);
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


}
