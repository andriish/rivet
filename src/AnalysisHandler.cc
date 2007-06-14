// -*- C++ -*-
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/RivetAIDA.hh"
#include "AIDA/IHistogramFactory.h"
#include "AIDA/IAnalysisFactory.h"
#include "AIDA/ITree.h"

using namespace std;


namespace Rivet {


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
    for (int i = 0, N = _analysisVector.size(); i < N; ++i) {
      _analysisVector[i]->finalize();
    }
    tree().commit();
  }


}
