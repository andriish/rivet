// -*- C++ -*-

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/RivetAIDA.hh"
//#include "AIDA/IHistogramFactory.h"
//#include "AIDA/IAnalysisFactory.h"
//#include "AIDA/ITree.h"
#include "AIDA/IManagedObject.h"
#include "AIDA/IDataPointSetFactory.h"


using namespace std;


namespace Rivet {


  inline Log& AnalysisHandler::getLog() { 
    return Log::getLog("Rivet.AnalysisHandler");
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
  }


  void AnalysisHandler::normalizeTree(AIDA::ITree& tree) {
    const vector<string> paths = tree.listObjectNames(".", true); // args set recursive listing
    getLog() << Log::INFO << "Number of objects in AIDA tree = " << paths.size() << endl;
    for (vector<string>::const_iterator path = paths.begin(); path != paths.end(); ++path) {
      AIDA::IManagedObject* obj = tree.find(*path);
      cout << "Testing type of " << *path << endl;
      if (obj) {
        AIDA::IHistogram1D* histo = (AIDA::IHistogram1D*) obj->cast("AIDA::IHistogram");
        if (histo) {
          cout << *path << " is a IH1D" << endl;
        }
      }
    }
  }


  AnalysisHandler::AnalysisHandler(string basefilename, HistoFormat storetype)
    : _nRun(0), _iRun(0) {
    _theAnalysisFactory = AIDA_createAnalysisFactory();
    setupFactories(basefilename, storetype);
  }


  AnalysisHandler::AnalysisHandler(AIDA::IAnalysisFactory& afac, string basefilename, HistoFormat storetype)
    : _nRun(0), _iRun(0), _theAnalysisFactory(&afac) {
    setupFactories(basefilename, storetype);
  }


  AnalysisHandler::~AnalysisHandler() {
    for (int i = 0, N = _analysisVector.size(); i < N; ++i) {
      delete _analysisVector[i]; 
    }
  }


  void AnalysisHandler::init(int i, int N) {
    _nRun = N;
    _iRun = i;
    for (int i = 0, N = _analysisVector.size(); i < N; ++i) {
      _analysisVector[i]->init();
      _analysisVector[i]->checkConsistency();
    }
  }


  void AnalysisHandler::analyze(const GenEvent & geneve) {
    Event event(geneve);
    for (int i = 0, N = _analysisVector.size(); i < N; ++i) {
      _analysisVector[i]->analyze(event);
    }
  }


  void AnalysisHandler::finalize() {
    Log& log = getLog();
    log << Log::INFO << "Finalising analysis" << endl;
    for (int i = 0, N = _analysisVector.size(); i < N; ++i) {
      _analysisVector[i]->finalize();
    }
    log << Log::INFO << "Normalising the AIDA tree" << endl;
    assert(_theTree != 0);
    normalizeTree(tree());
    tree().commit();
  }


}
