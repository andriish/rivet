// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TestMultiplicity class.
//

#include "TestMultiplicity.h"

using namespace Rivet;

TestMultiplicity::~TestMultiplicity() {}

void TestMultiplicity::init() {
  // Book histogram here.
}

void TestMultiplicity::analyze(const Event & event) {
  const FinalStateHCM & fs = event(fsproj);
  // int mult =
    fs.particles().size();
  // Fill histogram here.
}

void TestMultiplicity::finalize() {
  // Nothing to do here.
}


