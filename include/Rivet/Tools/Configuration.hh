// -*- C++ -*-
#ifndef RIVET_CONFIGURATION_HH 
#define RIVET_CONFIGURATION_HH 1

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/ParticleName.hh"


namespace Rivet {

  /// Typedef for a parameter name:value pair.
  typedef pair<string, string> Param;
  typedef const pair<string, string> cParam;

  /// Typedef for a parameter container.
  typedef map<string, string> ParamMap;
  typedef const map<string, string> cParamMap;

  /// @author Andy Buckley
  /// @date   2007-02-08
  class Configuration {
  public:
    /// Standard constructor
    Configuration() :
      numEvents(0), generatorName(""), beam1(Rivet::PROTON), beam2(Rivet::PROTON), 
      mom1(7000.0), mom2(7000.0), histoName("Rivet"), histoFormat(Rivet::AIDAML), 
      hepmlInFile(""), hepmlOutFile(""), hepmcInFile(""), hepmcOutFile(""), 
      useLogColors(true), runRivet(false), readHepMC(false), writeHepMC(false), 
      params(), analyses(), rngSeed(314159)
    { }


    /// Destructor
    ~Configuration() { };

  public:
    size_t numEvents;
    string generatorName;
    Rivet::ParticleName beam1, beam2;
    double mom1, mom2;
    string histoName;
    Rivet::HistoFormat histoFormat;
    string hepmlInFile, hepmlOutFile;
    string hepmcInFile, hepmcOutFile;
    Rivet::Log::LevelMap logLevels;
    bool useLogColors;
    bool runRivet;
    bool readHepMC, writeHepMC;
    ParamMap params;
    set<string> analyses;
    size_t rngSeed;
  };


  // Forward declarations
  class Analysis;
  namespace Commandline {
    Configuration parseArgs(size_t argc, char** argv);
  }


}

#endif
