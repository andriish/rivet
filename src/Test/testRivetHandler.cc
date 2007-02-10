// $Id: $

#include "Rivet/Rivet.hh"
#include "Rivet/RivetHandler.hh"
#include "Rivet/Analysis/Analysis.hh"

using namespace Rivet;

int main() {
  Analysis::getAnalysis(ANALYSIS_TEST).getInfo();

  RivetHandler rh;
  rh.addAnalysis(ANALYSIS_TEST);

  return EXIT_SUCCESS;
}
