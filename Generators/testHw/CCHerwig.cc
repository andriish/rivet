#include "Herwig_i/Herwig.h"
#include "HepMC/GenEvent.h"

int main() {
  Herwig foo("hwinstance");
  HepMC::GenEvent myevent;

  StatusCode status;

  status = foo.genInitialize();
  status = foo.callGenerator();
  status = foo.genFinalize();
  status = foo.fillEvt(&myevent);

  // Do something to store or display myevent using the HepMC interface...
  return EXIT_SUCCESS;
}
