// -*- C++ -*-

#include "Rivet/Analysis/RivetHandler.h"
#include "Rivet/Analysis/TestMultiplicity.h"

int main() {

  using namespace Rivet;

  RivetHandler rivet("test");
  rivet.addAnalysis(TestMultiplicity());

  rivet.init();

  rivet.finalize();

}

