#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include <limits>
#include <cmath>
#include <iostream>

using namespace std;


class Test : public Rivet::Analysis {
public:
  Test() : Analysis("Test") {}

  void init() {
    book(_h_test, "test", 50, 66.0, 116.0);
  }

  void analyze(const Rivet::Event & e) {
    cout << "Normal fill" << '\n';
    _h_test->fill(90., 1.);

    cout << "Underflow fill" << '\n';
    _h_test->fill(30.,1.);

    cout << "Overflow fill" << '\n';
    _h_test->fill(130.,1.);

     cout << "Inf fill" << '\n';
    try {
      _h_test->fill(numeric_limits<double>::infinity(), 1.);
    } catch (YODA::RangeError & e) {
      cerr << e.what() << '\n';
      if ( string(e.what()) != string("X is Inf") ) throw;
    }

    cout << "NaN fill" << '\n';
    try {
      _h_test->fill(numeric_limits<double>::quiet_NaN(), 1.);
    } catch (YODA::RangeError & e) {
      cerr << e.what() << '\n';
      if ( string(e.what()) != string("X is NaN") ) throw;
    }
  }

private:
  Rivet::Histo1DPtr _h_test;
};

DECLARE_RIVET_PLUGIN(Test);

int main(int argc, char* argv[]) {
  assert(argc > 1);

  Rivet::AnalysisHandler rivet;
  rivet.addAnalysis("Test");

  std::ifstream file(argv[1]);
  HepMC::IO_GenEvent hepmcio(file);
  HepMC::GenEvent* evt = hepmcio.read_next_event();
  double sum_of_weights = 0.0;
  if (! evt) {
  	cerr << "No events\n";
  	return 1;
  }

  while (evt) {
    // Analyse current event
    rivet.analyze(*evt);
    sum_of_weights += evt->weights()[0];

    // Clean up and get next event
    delete evt; evt = 0;
    hepmcio >> evt;
  }
  file.close();

  rivet.setCrossSection(1.0, 0.1);
  rivet.finalize();
  rivet.writeData("NaN.aida");

  return 0;
}
