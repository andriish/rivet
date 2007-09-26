// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/RivetAIDA.hh"
#include "AIDA/IManagedObject.h"

using namespace AIDA;


namespace Rivet {


  Log& AnalysisHandler::getLog() { 
    return Log::getLog("Rivet.AnalysisHandler");
  }


  AnalysisHandler& AnalysisHandler::addAnalysis(const string& analysisname) { 
    Analysis* analysis = AnalysisLoader::getAnalysis(analysisname);
    if (analysis) { // < Check for null analysis.
      analysis->_theHandler = this;
      _analyses.insert(analysis);
    }
    return *this;
  }


  void AnalysisHandler::setupFactories(string basefilename, HistoFormat storetype) {
    string filename(basefilename), storetypestr("");
    if (storetype == AIDAML) {
      filename += ".aida";
      storetypestr = "xml";
    } else if (storetype == FLAT) {
      filename += ".data";
      storetypestr = "flat";
    } else if (storetype == ROOT) {
      filename += ".root";
      storetypestr = "root";
    }
    _theTree = _theAnalysisFactory->createTreeFactory()->create(filename, storetypestr, false, true);
    _theHistogramFactory = _theAnalysisFactory->createHistogramFactory(tree());
    _theDataPointSetFactory = _theAnalysisFactory->createDataPointSetFactory(tree());
  }


  void AnalysisHandler::normalizeTree(ITree& tree) {
    const vector<string> paths = tree.listObjectNames("/", true); // args set recursive listing
    getLog() << Log::DEBUG << "Number of objects in AIDA tree = " << paths.size() << endl;
    const string tmpdir = "/RivetNormalizeTmp";
    tree.mkdir(tmpdir);
    for (vector<string>::const_iterator path = paths.begin(); path != paths.end(); ++path) {
      IManagedObject* obj = tree.find(*path);
      if (obj) {
        IHistogram1D* histo = dynamic_cast<IHistogram1D*>(obj);
        if (histo) {
          tree.mv(*path, tmpdir);
          const size_t lastslash = path->find_last_of("/");
          const string basename = path->substr(lastslash+1, path->length() - (lastslash+1));
          const string tmppath = tmpdir + "/" + basename;
          IHistogram1D* tmphisto = dynamic_cast<IHistogram1D*>(tree.find(tmppath));
          if (tmphisto) {
            getLog() << Log::DEBUG << "Temp histo " << tmppath << " exists" << endl;
            datapointsetFactory().create(*path, *tmphisto);
          }
          tree.rm(tmppath);
        }
      }
    }
    tree.rmdir(tmpdir);
  }


  AnalysisHandler::AnalysisHandler(string basefilename, HistoFormat storetype)
    : _nRun(0), _iRun(0) {
    _theAnalysisFactory = AIDA_createAnalysisFactory();
    setupFactories(basefilename, storetype);
  }


  AnalysisHandler::AnalysisHandler(IAnalysisFactory& afac, string basefilename, HistoFormat storetype)
    : _nRun(0), _iRun(0), _theAnalysisFactory(&afac) {
    setupFactories(basefilename, storetype);
  }


  void AnalysisHandler::init(int i, int N) {
    _nRun = N;
    _iRun = i;
    for (set<Analysis*>::iterator a = _analyses.begin(); a != _analyses.end(); ++a) {
      (*a)->init();
      (*a)->checkConsistency();
    }
  }


  void AnalysisHandler::analyze(const GenEvent & geneve) {
    Event event(geneve);
    for (set<Analysis*>::iterator a = _analyses.begin(); a != _analyses.end(); ++a) {
      (*a)->analyze(event);
    }
  }


  void AnalysisHandler::finalize() {
    Log& log = getLog();
    log << Log::INFO << "Finalising analysis" << endl;
    for (set<Analysis*>::iterator a = _analyses.begin(); a != _analyses.end(); ++a) {
      (*a)->finalize();
    }

    // Change AIDA histos into data point sets
    log << Log::INFO << "Normalising the AIDA tree" << endl;
    assert(_theTree != 0);
    normalizeTree(tree());
    tree().commit();

    // Delete analyses
    for (set<Analysis*>::iterator a = _analyses.begin(); a != _analyses.end(); ++a) {
      delete *a;
    }
    _analyses.clear();
    AnalysisLoader::closeAnalysisBuilders();
  }


}
