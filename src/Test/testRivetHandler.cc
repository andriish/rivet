// $Id: $

#include "Rivet/Rivet.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Analysis/Analysis.hh"

using namespace Rivet;

int main() {
  Analysis::getAnalysis(ANALYSIS_TEST);

  AnalysisHandler rh;
  rh.addAnalysis(ANALYSIS_TEST);

  return EXIT_SUCCESS;
}
