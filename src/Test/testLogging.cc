#include "Rivet/Logging.h"
#include <cstdlib>

using namespace std;
using namespace Rivet;

int main() {
  Logger& log = getLogger(); //("rivet.log");
  log.setPriority(log4cpp::Priority::INFO);
  log.info("This is some info");
	log.debug("This debug message will fail to write");
	log.alert("All hands abandon ship");
  
  return EXIT_SUCCESS;
}

