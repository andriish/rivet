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
    theTree = theAnalysisFactory->createTreeFactory()->create(filename, storetypestr, false, true);
    theHistogramFactory = theAnalysisFactory->createHistogramFactory(tree());
  }


  AnalysisHandler::AnalysisHandler(string basefilename, HistoFormat storetype)
    : nRun(0), iRun(0) {
    theAnalysisFactory = AIDA_createAnalysisFactory();
    setupFactories(basefilename, storetype);
  }


  AnalysisHandler::AnalysisHandler(AIDA::IAnalysisFactory& afac, string basefilename, HistoFormat storetype)
    : nRun(0), iRun(0), theAnalysisFactory(&afac) {
    setupFactories(basefilename, storetype);
  }


  AnalysisHandler::~AnalysisHandler() {
    for (int i = 0, N = analysisVector.size(); i < N; ++i) {
      delete analysisVector[i]; 
    }
  }


  void AnalysisHandler::init(int i, int N) {
    nRun = N;
    iRun = i;
    for (int i = 0, N = analysisVector.size(); i < N; ++i) {
      analysisVector[i]->init();
    }
  }


  void AnalysisHandler::analyze(const GenEvent & geneve) {
    Event event(geneve);
    for (int i = 0, N = analysisVector.size(); i < N; ++i) {
      analysisVector[i]->analyze(event);
    }
  }


  void AnalysisHandler::finalize() {
    for (int i = 0, N = analysisVector.size(); i < N; ++i) {
      analysisVector[i]->finalize();
    }
    tree().commit();
  }


//   RivetInfo AnalysisHandler::info() const {
//     RivetInfo ret;
//     for (int i = 0, N = analysisVector.size(); i < N; ++i) {
//       ret += analysisVector[i]->getInfo();
//     }
//     return ret;
//   }


}
