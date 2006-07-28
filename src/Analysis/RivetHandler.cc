// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RivetHandler class.
//

#include "Rivet/Analysis/RivetHandler.h"
#include "AIDA/ITreeFactory.h"

using namespace Rivet;

RivetHandler::RivetHandler() {}

RivetHandler::RivetHandler(string filename, string storetype,
			   AIDA::IAnalysisFactory & afac)
  : nRun(0), iRun(0), theAnalysisFactory(&afac) {
  theTree = afac.createTreeFactory()->create(filename, storetype, false, true);
  theHistogramFactory = afac.createHistogramFactory(tree());
}

RivetHandler::~RivetHandler() {
  for ( int i = 0, N = anaVector.size(); i < N; ++i ) delete anaVector[i];
}

void RivetHandler::init(int i, int N) {
  nRun = N;
  iRun = i;
  for ( int i = 0, N = anaVector.size(); i < N; ++i ) anaVector[i]->init();
}

void RivetHandler::analyze(const GenEvent & geneve) {
  Event event(geneve);
  for ( int i = 0, N = anaVector.size(); i < N; ++i )
    anaVector[i]->analyze(event);
}

void RivetHandler::finalize() {
  for ( int i = 0, N = anaVector.size(); i < N; ++i )
    anaVector[i]->finalize();
  tree().commit();
}

RivetInfo RivetHandler::info() const {
  RivetInfo ret;
  for ( int i = 0, N = anaVector.size(); i < N; ++i )
    ret += anaVector[i]->getInfo();
  return ret;
}
