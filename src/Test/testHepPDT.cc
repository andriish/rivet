#include "Rivet/Rivet.h"
#include "HepPDT/ParticleID.hh"
#include <cstdlib>

using namespace std;
using namespace Rivet;

int main() {
  Logger& log = getLogger();
  log.setPriority(LogPriority::INFO);
  
  HepPDT::ParticleID pInfo(11);
  log << LogPriority::INFO << "PID: " << pInfo.pid() << endlog;
  if (pInfo.isHadron()) log.info("It's a hadron");
  if (pInfo.isLepton()) log.info("It's a lepton");

  return EXIT_SUCCESS;
}

