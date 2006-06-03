// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AnalysisBase class.
//

#include "AnalysisBase.h"
#include "Rivet/Analysis/RivetHandler.h"

using namespace Rivet;

AnalysisBase::~AnalysisBase() {}

void AnalysisBase::init() {}

void AnalysisBase::analyze(const Event &) {}

void AnalysisBase::finalize() {}

RivetInfo AnalysisBase::getInfo() const {
  return info;
}

AIDA::IAnalysisFactory & AnalysisBase::analysisFactory() {
  return handler().analysisFactory();
}

AIDA::ITree & AnalysisBase::tree() {
  return handler().tree();
}

AIDA::IHistogramFactory & AnalysisBase::histogramFactory() {
  return handler().histogramFactory();
}

