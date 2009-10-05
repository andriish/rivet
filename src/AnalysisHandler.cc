// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/RivetAIDA.hh"
#include "LWH/AIManagedObject.h"

using namespace AIDA;

namespace Rivet {


  AnalysisHandler::AnalysisHandler(string basefilename, 
                                   string runname, HistoFormat storetype)
    : _runname(runname), _nRun(0), _iRun(0), _numEvents(0), _sumOfWeights(0.0) {
    _theAnalysisFactory = createAnalysisFactory();
    _setupFactories(basefilename, storetype);
  }


  AnalysisHandler::AnalysisHandler(IAnalysisFactory& afac, string basefilename, 
                                   string runname, HistoFormat storetype)
    : _runname(runname), _nRun(0), _iRun(0), _numEvents(0), _sumOfWeights(0.0), 
      _theAnalysisFactory(&afac) {
    _setupFactories(basefilename, storetype);
  }
  
  
  AnalysisHandler::~AnalysisHandler()
  {
  }


  Log& AnalysisHandler::getLog() { 
    return Log::getLog("Rivet.Analysis.Handler");
  }


  void AnalysisHandler::init(int i, int N) {
    getLog() << Log::DEBUG << "Initialising the analysis handler" << endl;
    _nRun = N;
    _iRun = i;
    _numEvents = 0;
    _sumOfWeights = 0.0;
    foreach (Analysis* a, _analyses) {
      getLog() << Log::DEBUG << "Initialising analysis: " << a->name() << endl;
      a->init();
      //getLog() << Log::DEBUG << "Checking consistency of analysis: " << a->name() << endl;
      //a->checkConsistency();
      getLog() << Log::DEBUG << "Done initialising analysis: " << a->name() << endl;
    }
    getLog() << Log::DEBUG << "Analysis handler initialised" << endl;
  }
  

  void AnalysisHandler::analyze(const GenEvent& ge) {
    Event event(ge);
    _numEvents++;
    _sumOfWeights += event.weight();
    #ifdef HEPMC_HAS_CROSS_SECTION
    if (ge.cross_section()) {
      /// @todo Use xs error?
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
    getLog() << Log::INFO << "Finalising analysis" << endl;
    foreach (Analysis* a, _analyses) {
      a->finalize();
    }

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


  /// Remove beam-incompatible analyses from the run list.
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
        IHistogram1D* histo = dynamic_cast<IHistogram1D*>(hobj);
        IProfile1D* prof = dynamic_cast<IProfile1D*>(hobj);
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
  
  
  std::vector<std::string> AnalysisHandler::analysisNames() {
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
    foreach (Analysis* a, _analyses) {
      a->setCrossSection(xs);
    }
    return *this;
  }

}
