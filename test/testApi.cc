#include "Rivet/AnalysisHandler.hh"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

using namespace std;

int main(int argc, char* argv[]) {
  assert(argc > 1);

  Rivet::AnalysisHandler ah;
  Rivet::Log::setLevel("Rivet", Rivet::Log::DEBUG);

  // Specify the analyses to be used
  ah.addAnalysis("EXAMPLE");
  ah.addAnalyses({{ "MC_JETS", "EXAMPLE_CUTS", "EXAMPLE_SMEAR" }});

  std::ifstream file(argv[1]);
  HepMC::IO_GenEvent hepmcio(file);
  HepMC::GenEvent* evt = hepmcio.read_next_event();
  if (! evt) {
  	cerr << "No events\n";
  	return 1;
  }

  while (evt) {
    // Analyse current event
    ah.analyze(*evt);

    // Clean up and get next event
    delete evt; evt = nullptr;
    hepmcio >> evt;
  }
  file.close();

  ah.setCrossSection(1.0, 0.1);

  ah.finalize();
  ah.writeData("out.yoda");

  return 0;
}
