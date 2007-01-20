// $Id: $


#include "LWH/AnalysisFactory.h"
#include "LWH/TreeFactory.h"
#include "LWH/Tree.h"
#include "LWH/HistogramFactory.h"
#include "LWH/Histogram1D.h"

#include <iostream>
#include <fstream>


using namespace LWH;
using namespace AIDA;

int main(int argc, char* argv[]) {
  IAnalysisFactory* af = new AnalysisFactory();
  ITree* tree =
    af->createTreeFactory()->create("test.data", "flat", false, true);
//   ITree* tree =
//     af->createTreeFactory()->create("test.aida", "xml", false, true);
  IHistogramFactory* hf = af->createHistogramFactory(*tree);
  IHistogram1D* h =
    hf->createHistogram1D("/test1", "Testing histo", 100, 0.0, 5.0);

  for (int i = 0 ; i < 10000 ; ++i) {
    double x = i/2000.0;
    double y = sin(2*x);
    //std::cout << "Filling: " << x << ", " << y << std::endl;
    h->fill(x, y);
  }

  tree->commit();

  return EXIT_SUCCESS;
}
