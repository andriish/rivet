// -*- C++ -*-
#ifndef RIVET_COMMANDLINE_HH 
#define RIVET_COMMANDLINE_HH 1

#include "Rivet/Rivet.hh"
#include "Rivet/HistoFormat.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Tools/Logging.hh"
#include "RivetGun/Configuration.hh"


namespace Rivet {
  namespace Commandline {

    Configuration parseArgs(size_t argc, char** argv);

  }  
}

#endif
