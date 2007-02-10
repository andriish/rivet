// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RivetHandler class.
//

#include "Rivet/RivetHandler.hh"
#include "Rivet/RivetAIDA.hh"
#include "AIDA/IHistogramFactory.h"
#include "AIDA/IAnalysisFactory.h"
#include "AIDA/ITree.h"

using namespace std;


namespace Rivet {

  void RivetHandler::setupFactories(string basefilename, HistoFormat storetype) {
    string filename(basefilename), storetypestr("");
    if (storetype == XML) {
      filename += ".aida";
      storetypestr = "xml";
    } else {
      filename += ".data";
      storetypestr = "flat";
    }
    theTree = theAnalysisFactory->createTreeFactory()->create(filename, storetypestr, false, true);
    theHistogramFactory = theAnalysisFactory->createHistogramFactory(*tree());
  }


  RivetHandler::RivetHandler(string basefilename, HistoFormat storetype)
    : nRun(0), iRun(0) {
    theAnalysisFactory = AIDA_createAnalysisFactory();
    setupFactories(basefilename, storetype);
  }


  RivetHandler::RivetHandler(AIDA::IAnalysisFactory& afac, string basefilename, HistoFormat storetype)
    : nRun(0), iRun(0), theAnalysisFactory(&afac) {
    setupFactories(basefilename, storetype);
  }


  RivetHandler::~RivetHandler() {
    /// @todo Can these for loops be made to use iterators?
    for ( int i = 0, N = analysisVector.size(); i < N; ++i ) {
      delete analysisVector[i]; 
    }
  }

  void RivetHandler::init(int i, int N) {
    nRun = N;
    iRun = i;
    /// @todo Can these for loops be made to use iterators?
    for ( int i = 0, N = analysisVector.size(); i < N; ++i ) {
      analysisVector[i]->init();
    }
  }

  void RivetHandler::analyze(const GenEvent & geneve) {
    Event event(geneve);
    /// @todo Can these for loops be made to use iterators?
    for ( int i = 0, N = analysisVector.size(); i < N; ++i ) {
      analysisVector[i]->analyze(event);
    }
  }

  void RivetHandler::finalize() {
    /// @todo Can these for loops be made to use iterators?
    for ( int i = 0, N = analysisVector.size(); i < N; ++i ) {
      analysisVector[i]->finalize();
    }
    if (tree()) { tree()->commit(); }
  }

  RivetInfo RivetHandler::info() const {
    RivetInfo ret;
    /// @todo Can these for loops be made to use iterators?
    for ( int i = 0, N = analysisVector.size(); i < N; ++i ) {
      ret += analysisVector[i]->getInfo();
    }
    return ret;
  }

}
