// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TestMultiplicity class.
//

#include "TestMultiplicity.h"

using namespace Rivet;
using namespace std;

TestMultiplicity::~TestMultiplicity() {}

void TestMultiplicity::init() {
  // Book histogram here.
}

void TestMultiplicity::analyze(const Event & event) {
  const FinalStateProjection & fs = event(fsproj);
  int mult =  fs.particles().size();
  cout << "Event multiplicity = " << mult << endl;
  // Fill histogram here.
}

void TestMultiplicity::finalize() {
  // Nothing to do here.
}


