// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AnalysisBase class.
//

#include "AnalysisBase.h"

using namespace Rivet;

AnalysisBase::~AnalysisBase() {}

void AnalysisBase::init() {}

void AnalysisBase::analyze(const Event & event) {}

void AnalysisBase::finalize() {}

RivetInfo AnalysisBase::getInfo() const {
  return info;
}

