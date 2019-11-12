#include <fstream>
#include "Rivet/Tools/RivetHepMC.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/AnalysisLoader.hh"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <hepmcfile> <ana1> [<ana2> ...]" << '\n';
    std::cout << "Available analyses:\n";
    for (const std::string& a : Rivet::AnalysisLoader::analysisNames())
      std::cout << "  " << a << "\n";
    std::cout << std::endl;
    return 1;
  }

  Rivet::AnalysisHandler ah;
  for (int i = 2; i < argc; ++i) {
    ah.addAnalysis(argv[i]);
  }

  std::shared_ptr<std::istream> istr;
  
  std::shared_ptr<Rivet::HepMC_IO_type> reader = Rivet::HepMCUtils::makeReader(argv[1], istr);
  
  std::shared_ptr<Rivet::RivetHepMC::GenEvent> evt = std::make_shared<Rivet::RivetHepMC::GenEvent>();

 
  while(reader && Rivet::HepMCUtils::readEvent(reader, evt)){
    ah.analyze(evt.get());
    evt.reset(new Rivet::RivetHepMC::GenEvent());
  }
  

  

  ah.setCrossSection(std::make_pair(1.0, 0.0));
  ah.finalize();
  ah.writeData("Rivet.yoda");

  return 0;
}
