// -*- C++ -*-
#include "RivetGun/RivetGun.hh"
#include "RivetGun/Configuration.hh"
#include "RivetGun/Commandline.hh"
#include "RivetGun/AvailableGenerators.hh"
#include "Rivet/Tools/Logging.hh"
using namespace Rivet;

#include "AGILe/AGILe.hh"
#include "AGILe/Particle.hh"
#include "AGILe/Generator.hh"

using namespace AGILe;

namespace Rivet {
  void generate(Generator* gen, Configuration& cfg, Log& log);
}

#ifdef FC_DUMMY_MAIN
extern "C" int F77_DUMMY_MAIN() { return 1; }
#endif


////////////////////////////////////////////////


int main(int argc, char* argv[]) {
  // Configuration flags
  Configuration cfg;
  try {
    cfg = Commandline::getConfigFromArgs(argc, argv, getAvailableGenNames());
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
    
    // Choose generator
    log << Log::INFO << "Requested generator = " << cfg.generatorName << endl;
    Generator* gen = makeNewGenerator(cfg.generatorName);
    log << Log::DEBUG << "Built new generator" << endl;
    if (!gen) {
      log << Log::ERROR << "No generator chosen... exiting" << endl;
      return EXIT_FAILURE;
    }
    
    // Configure and run the generator
    log << Log::DEBUG << "Start generation" << endl;
    generate(gen, cfg, log);
  }

  // Main loop exception handling.
  catch (runtime_error& e) { 
    cerr << "Error: " << e.what() << endl; 
    return EXIT_FAILURE;
  }


  return EXIT_SUCCESS;
}
