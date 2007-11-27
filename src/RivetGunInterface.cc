// -*- C++ -*-
#include "AGILe/AGILe.hh"
#include "AGILe/Particle.hh"
#include "AGILe/Generator.hh"
#include "AGILe/Loader.hh"
using namespace AGILe;

#include "RivetGun/RivetGun.hh"
#include "RivetGun/Configuration.hh"
#include "RivetGun/Commandline.hh"
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
  
  // Configuration flags
  Configuration cfg;
  try {
    // Use the list of available generators to configure 
    // and use the command line handler.
    cfg = Commandline::getConfigFromArgs(argc, argv, Loader::getAvailableGens());
  } catch (TCLAP::ArgException& e) { 
    cerr << "Command line error: " << e.error() << " for arg " << e.argId() << endl; 
    return EXIT_FAILURE;
  } catch (runtime_error& e) { 
    cerr << "Error: " << e.what() << endl; 
    return EXIT_FAILURE;
  }


  try {
    // Set log levels from command line and get a logger
    Log::setDefaultLevels(cfg.logLevels);
    Log& log = Log::getLog("RivetGun.Main");


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
