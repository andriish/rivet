// -*- C++ -*-
#ifndef RIVETGUN_COMMANDLINE_HH 
#define RIVETGUN_COMMANDLINE_HH 1

#include "RivetGun/RivetGun.hh"
#include "RivetGun/Configuration.hh"
#include <tclap/CmdLine.h>


namespace Rivet {

  namespace Commandline {
    /// @author Andy Buckley
    /// @date   2007-02-08


    /// Parse the command line arguments and make a Configuration object
    Configuration getConfigFromArgs(int argc, char* argv[], vector<string> genNames);


    /// Get param defn lines from an istream and add them into the supplied map
    void handleConfigStream(std::istream& in, map<string, string>& pmap);


    /// Add command line arguments to specify the generator to be used.
    void addGeneratorArgs(TCLAP::CmdLine& cmd, vector<string>& genNames);


    /// Add command line parsing of the beam and momenta parameters.
    void addInitialStateArgs(TCLAP::CmdLine& cmd);


    /// Add command line argument handling for generator parameters.
    void addParamArgs(TCLAP::CmdLine& cmd);


    /// Set the generator components of the Configuration object based on the
    /// command line arguments.  Be careful - this function deletes the argument
    /// objects.
    void useGeneratorArgs(TCLAP::CmdLine& cmd, Configuration& config);


    /// Set the initial state (beam and momenta) components of the Configuration
    /// object based on the command line arguments.  Be careful - this function
    /// deletes the argument objects.
    void useInitialStateArgs(TCLAP::CmdLine& cmd, Configuration& config);

    
    /// Set the generator parameter components of the Configuration object based
    /// on the command line arguments.  Be careful - this function deletes the
    /// argument objects.
    void useParamArgs(TCLAP::CmdLine& cmd, Configuration& config);
    

  }  
}

#endif
