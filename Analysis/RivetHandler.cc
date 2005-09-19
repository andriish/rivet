// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RivetHandler class.
//

#include "RivetHandler.h"

using namespace Rivet;

RivetHandler::~RivetHandler() {
  for ( int i = 0, N = anaVector.size(); i < N; ++i ) delete anaVector[i];
}

void RivetHandler::init(int i, int N) {
  nRun = N;
  iRun = i;
  // *** ATTENTION *** here we should initialize the histogram factory.
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
  // *** ATTENTION *** here we should write out the histograms from
  // the histogram factory.
}

