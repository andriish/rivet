// -*- C++ -*-
#include "AGILe/AGILe.hh"
//#include "AGILe/Particle.hh"
//#include "AGILe/Generator.hh"
//using AGILe::Generator;
//using AGILe::Loader;

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Configuration.hh"
#include "Rivet/Tools/Logging.hh"
using namespace Rivet;

#ifdef FC_DUMMY_MAIN
extern "C" int F77_DUMMY_MAIN() { return 1; }
#endif


////////////////////////////////////////////////


// Forward declare the function that actually runs 
// the event loop and analysis, which is isolated from 
// the administration and param parsing code in main().
namespace Rivet {
  void generate(Configuration& cfg, Log& log);
}


// The main function of the rivetgun executable.
int main(int argc, char* argv[]) {
  
  // Parse command line into a configuration object.
  Configuration cfg;
  try {
    cfg = Commandline::parseArgs(argc, argv);
  } catch (runtime_error& e) { 
    cerr << "Error: " << e.what() << endl; 
    return EXIT_FAILURE;
  }


  // Set log levels from command line and get a logger
  Log::setDefaultLevels(cfg.logLevels);
  Log::setUseColors(cfg.useLogColors);
  //AGILe::Log::setDefaultLevels(cfg.logLevels);
  AGILe::Log::setUseColors(cfg.useLogColors);
  Log& log = Log::getLog("RivetGun.Main");


  try {
    // Configure and run the generator.
    generate(cfg, log);
  } 
  
  // Main loop exception handling.
  catch (runtime_error& e) { 
    cerr << "Error: " << e.what() << endl; 
    return EXIT_FAILURE;
  }


  return EXIT_SUCCESS;
}
