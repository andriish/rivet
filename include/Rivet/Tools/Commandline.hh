// -*- C++ -*-
#ifndef RIVET_COMMANDLINE_HH 
#define RIVET_COMMANDLINE_HH 1

#include "Rivet/Rivet.hh"
#include "Rivet/HistoFormat.hh"
#include "Rivet/Analysis/AnalysisLoader.hh"
#include "Rivet/Tools/Logging.hh"
#include <tclap/CmdLine.h>


namespace Rivet {

  class Analysis;


  namespace Commandline {
    /// @author Andy Buckley
    /// @date   2007-02-08


    /// Add command line parsing of analyses and their negations.
    void addAnalysisArgs(TCLAP::CmdLine& cmd,
                         TCLAP::MultiArg<string>*& analysesArg,
                         TCLAP::SwitchArg*& analysesAllArg,
                         TCLAP::ValuesConstraint<string>*& anaNameConstraint);

    /// Add command line parsing of histogram file and format specifications.
    void addHistoArgs(TCLAP::CmdLine& cmd,
                      TCLAP::ValueArg<string>*& histoNameArg,
                      TCLAP::ValueArg<string>*& histoTypeArg,
                      TCLAP::ValuesConstraint<string>*& histoTypeConstraint);

    /// Add command line options to control the logging system.
    void addLoggingArgs(TCLAP::CmdLine& cmd, 
                        TCLAP::MultiArg<string>*& logsArg);
    
    /// Set the Rivet analysis components of the Configuration object based on
    /// the command line arguments.  Be careful - this function deletes the
    /// argument objects.
    void useAnalysisArgs(TCLAP::CmdLine& cmd,
                         TCLAP::MultiArg<string>* analysesArg,
                         TCLAP::SwitchArg* analysesAllArg,
                         set<string>& cfgAnalyses);

    /// Set the histogramming components of the Configuration object based on
    /// the command line arguments.  Be careful - this function deletes the
    /// argument objects.
    void useHistoArgs(TCLAP::CmdLine& cmd,
                      TCLAP::ValueArg<string>* histoNameArg,
                      TCLAP::ValueArg<string>* histoTypeArg,
                      string& cfgHistoFileName,
                      HistoFormat& cfgHistoFormat);

    /// Set the logging components of the Configuration object based on the
    /// command line arguments.  Be careful - this function deletes the argument
    /// object.
    void useLoggingArgs(TCLAP::CmdLine& cmd,
                        TCLAP::MultiArg<string>* logsArg,
                        Log::LevelMap& cfgLogLevels);

  }  
}

#endif
