// -*- C++ -*-
#include "AGILe/AGILe.hh"
#include "AGILe/Particle.hh"
#include "AGILe/Generator.hh"
#include "AGILe/Loader.hh"
using namespace AGILe;

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
  void generate(Generator* gen, Configuration& cfg, Log& log);
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
  Log& log = Log::getLog("RivetGun.Main");


  try {
    // Load generator libraries
    Loader::initialize();
    log << Log::INFO << "Requested generator = " << cfg.generatorName << endl;
    try {
      Loader::loadGenLibs(cfg.generatorName);
    } catch (runtime_error& e) {
      log << Log::ERROR << "Error when loading the generator library:" 
          << endl << e.what() << endl;
    }

    // Make a generator object.
    Generator* gen = Loader::createGen();
    if (!gen) {
      log << Log::ERROR << "No generator chosen... exiting" << endl;
      return EXIT_FAILURE;
    }

    // Configure and run the generator.
    generate(gen, cfg, log);

    // Shut down dynamic loader.
    Loader::destroyGen(gen);
    Loader::finalize();  
  } 
  
  // Main loop exception handling.
  catch (runtime_error& e) { 
    cerr << "Error: " << e.what() << endl; 
    return EXIT_FAILURE;
  }


  return EXIT_SUCCESS;
}
