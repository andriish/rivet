#include "Rivet/AnalysisHandler.hh"
#include "Rivet/AnalysisLoader.hh"
//#include "HepMC/IO_GenEvent.h"
#include "HepMCContrib/ReaderFactory.h"

using namespace std;

int main(int argc, char** argv) {
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " <hepmcfile> <ana1> [<ana2> ...]" << endl;
    return 1;
  }

  foreach (const string& a, Rivet::AnalysisLoader::analysisNames())
    cout << a << endl;

  Rivet::AnalysisHandler ah;
  for (int i = 2; i < argc; ++i) {
    ah.addAnalysis(argv[i]);
  }

  std::shared_ptr<HepMC::Reader> reader = HepMC::make_reader(argv[1]);
  
  while(!reader->failed()){
    
    HepMC::GenEvent evt;
    reader->read_event(evt);
    ah.analyze(evt);
  }
  
  reader->close();

  ah.setCrossSection(1.0);
  ah.finalize();
  ah.writeData("Rivet.yoda");

  return 0;
}
