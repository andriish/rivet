#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "HepPDT/ParticleID.hh"
#include <cstdlib>

using namespace std;
using namespace Rivet;

int main() {
  Log& log = Log::getLog("Test");
  
  HepPDT::ParticleID pInfo(11);
  log << Log::INFO << "PID: " << pInfo.pid() << endl;
  if (pInfo.isHadron()) log.info("It's a hadron");
  if (pInfo.isLepton()) log.info("It's a lepton");

  return EXIT_SUCCESS;
}

