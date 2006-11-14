#include "Rivet/Logging.h"
#include <cstdlib>

using namespace std;
using namespace Rivet;

int main() {
  Logger& log = getLogger();
  log.setPriority(LogPriority::INFO);
  log.info("This is some info");
	log.debug("This debug message will fail to write");
	log.alert("All hands abandon ship");
  
  // Check that we don't try to instantiate the logger twice
  Logger& log2 = getLogger();
  log2 << LogPriority::ALERT << "Another alert!" << endlog;

  return EXIT_SUCCESS;
}

