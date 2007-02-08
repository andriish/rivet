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


namespace Rivet {

  /// @todo Abstract histo booking in a function elsewhere. Auto-append file extension.

  RivetHandler::RivetHandler() {
    theAnalysisFactory = AIDA_createAnalysisFactory();
    theTree = theAnalysisFactory->createTreeFactory()->create("Rivet.data", "flat", false, true);
    //ITree* tree = theAnalysisFactory->createTreeFactory()->create("Rivet.aida.xml", "xml", false, true);
    theHistogramFactory = theAnalysisFactory->createHistogramFactory(*theTree);
  }

  RivetHandler::RivetHandler(string filename, string storetype, AIDA::IAnalysisFactory& afac)
    : nRun(0), iRun(0), theAnalysisFactory(&afac) {
    theTree = afac.createTreeFactory()->create(filename, storetype, false, true);
    theHistogramFactory = afac.createHistogramFactory(*tree());
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
