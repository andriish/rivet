// $Id: $

#include "Rivet/RivetAIDA.hh"
#include "AIDA/IHistogram1D.h"
using namespace AIDA;

#include <iostream>
#include <fstream>


/// Plugin function to return an AIDA system (LWH impl.)
#include "LWH/AnalysisFactory.h"
extern "C" AIDA::IAnalysisFactory* AIDA_createAnalysisFactory() {
  return new LWH::AnalysisFactory();
}


/// Test the histogramming
int main(int argc, char* argv[]) {
  IAnalysisFactory* af = AIDA_createAnalysisFactory();
  //ITree* tree = af->createTreeFactory()->create("test.data", "flat", false, true);
  ITree* tree = af->createTreeFactory()->create("test.aida.xml", "xml", false, true);
  IHistogramFactory* hf = af->createHistogramFactory(*tree);
  IHistogram1D* h1 = hf->createHistogram1D("/test1", "Testing histo", 100, 0.0, 5.0);

  for (int i = 0 ; i < 10000 ; ++i) {
    double x = i/2000.0;
    double y = sin(x);
    //std::cout << "Filling: " << x << ", " << y << std::endl;
    h1->fill(x, y*y);
  }

  tree->commit();

  return EXIT_SUCCESS;
}
