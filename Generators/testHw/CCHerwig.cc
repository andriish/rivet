#include "Herwig_i/Herwig.h"
#include "HepMC/GenEvent.h"

struct
{
  int istg;
} evtcon_;


int main() {
  Herwig foo("hwinstance");
  HepMC::GenEvent myevent;

  evtcon_.istg = 0;
  StatusCode status;

  status = foo.genInitialize();
  status = foo.callGenerator();
  status = foo.genFinalize();
  status = foo.fillEvt(&myevent);

  // Do something to store or display myevent using the HepMC interface...
  return EXIT_SUCCESS;
}
