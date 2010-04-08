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


  AnalysisHandler::AnalysisHandler(string basefilename,
                                   string runname, HistoFormat storetype)
    : _runname(runname), _numEvents(0), 
      _sumOfWeights(0.0), _xs(-1.0), _initialised(false)
  {
    _theAnalysisFactory = createAnalysisFactory();
    _setupFactories(basefilename, storetype);
    initializeParticleNames();
  }


  AnalysisHandler::AnalysisHandler(IAnalysisFactory& afac, string basefilename,
                                   string runname, HistoFormat storetype)
    : _runname(runname), _numEvents(0), 
      _sumOfWeights(0.0), _xs(-1.0), _initialised(false), 
      _theAnalysisFactory(&afac) 
  {
    _setupFactories(basefilename, storetype);
    initializeParticleNames();
  }


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
    foreach (const Analysis* a, analyses()) {
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
    
    foreach (Analysis* a, _analyses) {
      getLog() << Log::DEBUG << "Initialising analysis: " << a->name() << endl;
      // Allow projection registration in the init phase onwards
      a->_allowProjReg = true;
      a->init();
      //getLog() << Log::DEBUG << "Checking consistency of analysis: " << a->name() << endl;
      //a->checkConsistency();
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
//     // Ensure that beam details match those from first event @todo This is too expensive
//     const BeamPair beams = Rivet::beamIds(ge);
//     const double sqrts = Rivet::sqrtS(ge);
//     if (!compatible(beams, _beams) || !fuzzyEquals(sqrts, sqrtS())) {
//       getLog() << Log::ERROR << "Event beams mismatch: "
//                << beams << " @ " << sqrts/GeV << " GeV" << " vs. first beams "
//                << this->beams() << " @ " << this->sqrtS()/GeV << " GeV" << endl;
//       exit(1);
//     }

    
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
    foreach (Analysis* a, _analyses) {
      getLog() << Log::DEBUG << "About to run analysis " << a->name() << endl;
      a->analyze(event);
      getLog() << Log::DEBUG << "Finished running analysis " << a->name() << endl;
    }
  }


  void AnalysisHandler::finalize() {
    assert(_initialised);
    getLog() << Log::INFO << "Finalising analyses" << endl;
    foreach (Analysis* a, _analyses) {
      a->finalize();
    }

    // Print out number of events processed
    getLog() << Log::INFO << "Processed " << _numEvents << " event" << (_numEvents == 1 ? "" : "s") << endl;

    // Change AIDA histos into data point sets
    getLog() << Log::DEBUG << "Converting histograms to scatter plots" << endl;
    assert(_theTree != 0);
    _normalizeTree(tree());

    // Delete analyses
    getLog() << Log::DEBUG << "Deleting analyses" << endl;
    foreach (Analysis* a, _analyses) {
      delete a;
    }
    _analyses.clear();

    // Delete singletons
    //ProjectionHandler::destroy();
    //HistoHandler::destroy();

    // Print out MCnet boilerplate
    cout << endl;
    cout << "The MCnet usage guidelines apply to Rivet: see http://www.montecarlonet.org/GUIDELINES" << endl;
    cout << "Please acknowledge plots made with Rivet analyses, and cite arXiv:1003.0694 (http://arxiv.org/abs/1003.0694)" << endl;
  }


  AnalysisHandler& AnalysisHandler::addAnalysis(const string& analysisname) {
    Analysis* analysis = AnalysisLoader::getAnalysis(analysisname);
    if (analysis) { // < Check for null analysis.
      getLog() << Log::DEBUG << "Adding analysis '" << analysisname << "'" << endl;
      analysis->_analysishandler = this;
      _analyses.insert(analysis);
    }
    return *this;
  }


  AnalysisHandler& AnalysisHandler::removeAnalysis(const string& analysisname) {
    Analysis* toremove = 0;
    foreach (Analysis* a, _analyses) {
      if (a->name() == analysisname) {
        toremove = a;
        break;
      }
    }
    if (toremove) {
      getLog() << Log::DEBUG << "Removing analysis '" << analysisname << "'" << endl;
      _analyses.erase(toremove);
      delete toremove;
    }
    return *this;
  }


  AnalysisHandler& AnalysisHandler::removeIncompatibleAnalyses(const BeamPair& beams) {
    vector<Analysis*> todelete;
    foreach (Analysis* a, _analyses) {
      if (! a->isCompatible(beams)) {
        todelete.push_back(a);
      }
    }
    foreach (Analysis* a, todelete) {
      getLog() << Log::WARN << "Removing incompatible analysis '"
               << a->name() << "'" << endl;
      _analyses.erase(a);
      delete a;
    }
    return *this;
  }


  void AnalysisHandler::_setupFactories(string basefilename, HistoFormat storetype) {
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
    _theTree = _theAnalysisFactory->createTreeFactory()->create(filename, storetypestr, false, true);
    _theHistogramFactory = _theAnalysisFactory->createHistogramFactory(tree());
    _theDataPointSetFactory = _theAnalysisFactory->createDataPointSetFactory(tree());
  }


  void AnalysisHandler::commitData() {
    tree().commit();
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
    foreach (Analysis* a, _analyses) {
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
    foreach (const Analysis* a, _analyses) {
      if (!rtn) rtn = a->needsCrossSection();
      if (rtn) break;
    }
    return rtn;
  }


  AnalysisHandler& AnalysisHandler::setCrossSection(double xs) {
    _xs = xs;
    foreach (Analysis* a, _analyses) {
      a->setCrossSection(xs);
    }
    return *this;
  }


  bool AnalysisHandler::hasCrossSection() const {
    return (crossSection() >= 0);
  }

  AnalysisHandler& AnalysisHandler::addAnalysis(Analysis* analysis) {
    analysis->_analysishandler = this;
    _analyses.insert(analysis);
    return *this;
  }

  BeamPair AnalysisHandler::beamIds() const { 
    return Rivet::beamIds(beams());
  }

  double AnalysisHandler::sqrtS() const {
    return Rivet::sqrtS(beams());
  }


  void AnalysisHandler::initializeParticleNames() {
    Rivet::s_pnames[ELECTRON] = "ELECTRON";
    Rivet::s_pnames[POSITRON] = "POSITRON";
    Rivet::s_pnames[PROTON] = "PROTON";
    Rivet::s_pnames[ANTIPROTON] = "ANTIPROTON";
    Rivet::s_pnames[PHOTON] = "PHOTON";
    Rivet::s_pnames[NEUTRON] = "NEUTRON";
    Rivet::s_pnames[ANTINEUTRON] = "ANTINEUTRON";
    Rivet::s_pnames[MUON] = "MUON";
    Rivet::s_pnames[ANTIMUON] = "ANTIMUON";
    Rivet::s_pnames[NU_E] = "NU_E";
    Rivet::s_pnames[NU_EBAR] = "NU_EBAR";
    Rivet::s_pnames[NU_MU] = "NU_MU";
    Rivet::s_pnames[NU_MUBAR] = "NU_MUBAR";
    Rivet::s_pnames[NU_TAU] = "NU_TAU";
    Rivet::s_pnames[NU_TAUBAR] = "NU_TAUBAR";
    Rivet::s_pnames[PIPLUS] = "PIPLUS";
    Rivet::s_pnames[PIMINUS] = "PIMINUS";
    Rivet::s_pnames[TAU] = "TAU";
    Rivet::s_pnames[WPLUSBOSON] = "WPLUSBOSON";
    Rivet::s_pnames[WMINUSBOSON] = "WMINUSBOSON";
    Rivet::s_pnames[ZBOSON] = "ZBOSON";
    Rivet::s_pnames[HIGGS] = "HIGGS";
    Rivet::s_pnames[ANTITAU] = "ANTITAU";
    Rivet::s_pnames[PHOTOELECTRON] = "PHOTOELECTRON";
    Rivet::s_pnames[PHOTOPOSITRON] = "PHOTOPOSITRON";
    Rivet::s_pnames[PHOTOMUON] = "PHOTOMUON";
    Rivet::s_pnames[PHOTOANTIMUON] = "PHOTOANTIMUON";
    Rivet::s_pnames[PHOTOTAU] = "PHOTOTAU";
    Rivet::s_pnames[PHOTOANTITAU] = "PHOTOANTITAU";
    Rivet::s_pnames[ANY] = "*";
  }



}
