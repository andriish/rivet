#include "Rivet/Logging.h"
#include <cstdlib>

using namespace std;
using namespace Rivet;

int main() {
  Logger& log = getLogger(); //("rivet.log");
  log.setPriority(LogPriority::INFO);
  log.info("This is some info");
	log.debug("This debug message will fail to write");
	log.alert("All hands abandon ship");
  
  log << LogPriority::ALERT << "Another alert!" << endlog;

  return EXIT_SUCCESS;
}

