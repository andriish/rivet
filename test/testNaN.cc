#include "Rivet/AnalysisHandler.hh"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "Rivet/Analysis.hh"
#include "Rivet/RivetYODA.hh"

using namespace std;

class Test : public Rivet::Analysis {
public:
  Test() : Analysis("Test") {}

  void init() {
    _h_test = bookHisto1D("test", 50, 66.0, 116.0);
  }

  void analyze(const Rivet::Event & e) {
    _h_test->fill(90.,1.);
    _h_test->fill(1./0.,1.);
    _h_test->fill(sqrt(-1.),1.);
    _h_test->fill(30.,1.);
    _h_test->fill(130.,1.);
  }

private:
  Rivet::Histo1DPtr _h_test;
};

DECLARE_RIVET_PLUGIN(Test);

int main() {
  Rivet::AnalysisHandler rivet;
  rivet.addAnalysis("Test");

  std::istream* file = new std::fstream("testApi.hepmc", std::ios::in);
  HepMC::IO_GenEvent hepmcio(*file);
  HepMC::GenEvent* evt = hepmcio.read_next_event();
  double sum_of_weights = 0.0;
  while (evt) {
    // Analyse current event
    rivet.analyze(*evt);
    sum_of_weights += evt->weights()[0];

    // Clean up and get next event
    delete evt; evt = 0;
    hepmcio >> evt;
  }
  delete file; file = 0;

  rivet.setCrossSection(1.0);
  rivet.finalize();
  rivet.writeData("NaN.aida");

  return 0;
}
