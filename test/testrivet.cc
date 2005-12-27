// -*- C++ -*-

#include "Rivet/Analysis/RivetHandler.h"
#include "Rivet/Analysis/Examples/TestMultiplicity.h"

int main() {

  using namespace Rivet;

  RivetHandler rivet("test");
  rivet.addAnalysis(TestMultiplicity());

  rivet.init();

  cout << rivet.info();

  rivet.finalize();

}

