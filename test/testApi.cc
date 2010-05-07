#include "Rivet/AnalysisHandler.hh"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

using namespace std;

int main() {
  // Old type of constructor
  /// @deprecated Use new kind which specifies output names only when writing
  Rivet::AnalysisHandler rivet_old("out", "", Rivet::AIDAML);

  // New type
  Rivet::AnalysisHandler rivet;

  // Specify the analyses to be used
  rivet.addAnalysis("D0_2008_S7554427");
  vector<string> moreanalyses(1, "D0_2007_S7075677");
  rivet.addAnalyses(moreanalyses);

  rivet.init(); //< Obsolete, but allowed for compatibility

  std::istream* file = new std::fstream("testApi.hepmc", std::ios::in);
  HepMC::IO_GenEvent hepmcio(*file);
  HepMC::GenEvent* evt = hepmcio.read_next_event();
  double sum_of_weights = 0.0;
  while (evt) {
    rivet.analyze(*evt);
    // Problem with HepMC file portability: temporarily disable
    // sum_of_weights += evt->weights()[0];

    // Clean up and get next event
    delete evt; evt = 0;
    hepmcio >> evt;
  }
  delete file; file = 0;

  rivet.setCrossSection(1.0);
  rivet.setSumOfWeights(sum_of_weights); // not necessary, but allowed
  rivet.finalize();
  rivet.writeData("out");

  return 0;
}
