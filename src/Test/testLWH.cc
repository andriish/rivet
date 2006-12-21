// $Id: $


#include "LWH/AnalysisFactory.h"
#include "LWH/HistogramFactory.h"
#include "LWH/Histogram1D.h"

#include <iostream>

using namespace LWH;

int main(int argc, char* argv[]) {
  //IAnalysisFactory* af = new AnalysisFactory();
  //IHistogramFactory* hf = af->createHistogramFactory();
  //IHistogram1D* h = hf->createHistogram1D("test1", "Testing histo", 10, 0.0, 10.0);

  Histogram1D h(20, 0.0, 5.0);

  for (int i = 0 ; i < 10000 ; ++i) {
    h.fill(i/2000.0);
  }

  h.writeXML(std::cout, "/test1", "foooo");

  return EXIT_SUCCESS;
}
