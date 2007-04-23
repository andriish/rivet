#include "Rivet/Tools/Logging.hh"
#include <cstdlib>

using namespace std;
using namespace Rivet;

int main() {
  Log::LevelMap levels;
  levels["Rivet.Test"] = Log::DEBUG;
  levels["Rivet.Test.Foo"] = Log::INFO;
  Log::setDefaultLevels(levels);

  Log& log1 = Log::getLog("Rivet.Test.Foo");
  Log& log2 = Log::getLog("Rivet.Test.Bar");

  log1.info("This is some info");
	log1.debug("This debug message will fail to write");
	log1.warn("All hands abandon ship");

  log2.info("This is some more info");
	log2.debug("This debug message should write out");
	log2.warn("Another warning...");
  
  return EXIT_SUCCESS;
}

