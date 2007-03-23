#include "Rivet/Tools/Logging.hh"
#include <cstdlib>

using namespace std;
using namespace Rivet;

int main() {
  Log& log = Log::getLog("RivetTest");
  log.info("This is some info");
	log.debug("This debug message will fail to write");
	log.warn("All hands abandon ship");
  
  /// @todo Check that we don't try to instantiate the logger twice
  Log& log2 = Log::getLog("RivetTest2");
  log2 << Log::WARN << "Another alert!" << endl;

  return EXIT_SUCCESS;
}

