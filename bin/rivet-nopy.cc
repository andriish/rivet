#include "Rivet/AnalysisHandler.hh"
#include "Rivet/AnalysisLoader.hh"
//#include "HepMC/IO_GenEvent.h"
#include "Rivet/Tools/RivetHepMC.hh"

using namespace std;

int main(int argc, char** argv) {
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " <hepmcfile> <ana1> [<ana2> ...]" << '\n';
    cout << "Available analyses:\n";
    for (const string& a : Rivet::AnalysisLoader::analysisNames())
      cout << "  " << a << "\n";
    cout << endl;
    return 1;
  }

  Rivet::AnalysisHandler ah;
  for (int i = 2; i < argc; ++i) {
    ah.addAnalysis(argv[i]);
  }

  std::ifstream istr(argv[1], std::ios::in);
  
  std::shared_ptr<Rivet::HepMC_IO_type> reader = Rivet::HepMCUtils::makeReader(istr);
  
  std::shared_ptr<Rivet::RivetHepMC::GenEvent> evt = make_shared<Rivet::RivetHepMC::GenEvent>();
  
  while(reader && Rivet::HepMCUtils::readEvent(reader, evt)){
    ah.analyze(evt.get());
    evt.reset(new Rivet::RivetHepMC::GenEvent());
  }
  
  istr.close();

  ah.setCrossSection(1.0, 0.0);
  ah.finalize();
  ah.writeData("Rivet.yoda");

  return 0;
}
