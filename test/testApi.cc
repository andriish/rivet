#include "Rivet/AnalysisHandler.hh"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

using namespace std;

int main() {
  string fname = "out";
  Rivet::AnalysisHandler rivet(fname, "", Rivet::AIDAML);

  // specify the analyses to be used
  rivet.addAnalysis("D0_2008_S7554427");
  vector<string> moreanalyses(1, "D0_2007_S7075677");
  rivet.addAnalyses(moreanalyses);

  rivet.init(); // obsolete, but allowed for compatibility

  std::istream* file = new std::fstream("testApi.hepmc", std::ios::in);
  HepMC::IO_GenEvent hepmcio(*file);
  HepMC::GenEvent* evt = hepmcio.read_next_event();
  double sum_of_weights = 0.0;
  while (evt) {
    rivet.analyze(*evt);
    sum_of_weights+=evt->weights()[0];
    // clean up and get next event
    delete evt;
    hepmcio >> evt;
  }
  delete file;

  rivet.setCrossSection(1.0);
  rivet.setSumOfWeights(sum_of_weights); // not necessary, but allowed
  rivet.finalize();
  rivet.commitData();

  return 0;
}
