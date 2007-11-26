// -*- C++ -*-
#include "RivetGun/RivetGun.hh"
#include "RivetGun/Commandline.hh"
#include "RivetGun/Configuration.hh"
#include "Rivet/AnalysisHandler.fhh"
#include "Rivet/Tools/Commandline.hh"
#include <tclap/CmdLine.h>
#include <fstream>


namespace Rivet {

  class Analysis;

  namespace Commandline {

    // Args as static variables
    static TCLAP::ValueArg<string> *genNameArg(0), *genFileArg(0);
    static TCLAP::ValuesConstraint<string>* genNameConstraint(0);
    static TCLAP::ValueArg<size_t>* rngSeedArg(0);
    static TCLAP::ValueArg<string> *beam1Arg(0), *beam2Arg(0);
    static TCLAP::ValueArg<double> *mom1Arg(0), *mom2Arg(0);
    static TCLAP::ValuesConstraint<string>* beamNameConstraint(0);
    static TCLAP::MultiArg<string>* paramsArg(0);
    static TCLAP::MultiArg<string>* paramsFilesArg(0);
    static TCLAP::MultiArg<string>* analysesArg(0);
    static TCLAP::SwitchArg* analysesAllArg(0);
    static TCLAP::ValuesConstraint<string>* anaNameConstraint(0);
    static TCLAP::ValueArg<string>* histoNameArg(0);
    static TCLAP::ValueArg<string>* histoTypeArg(0);
    static TCLAP::ValuesConstraint<string>* histoTypeConstraint(0);    
    static TCLAP::MultiArg<string>* logsArg(0);


    /// Clear all the arg object pointers
    void discardArgObjects() {
      if (genNameArg) { delete genNameArg; genNameArg = 0; }
      if (genFileArg) { delete genFileArg; genFileArg = 0; }
      if (rngSeedArg) { delete rngSeedArg; rngSeedArg = 0; }
      if (beam1Arg) { delete beam1Arg; beam1Arg = 0; }
      if (mom1Arg) { delete mom1Arg; mom1Arg = 0; }
      if (beam2Arg) { delete beam2Arg; beam2Arg = 0; }
      if (mom2Arg) { delete mom2Arg; mom2Arg = 0; }
      if (paramsArg) { delete paramsArg; paramsArg = 0; }
      if (paramsFilesArg) { delete paramsFilesArg; paramsFilesArg = 0; }
    }


    // Build a Configuration object from the command line arguments, using TCLAP
    Configuration getConfigFromArgs(int argc, char* argv[], vector<string> genNames) {
      TCLAP::CmdLine cmd("Runs an event generator using the RivetGun and AGILe interface libraries", ' ', "1.0"); 

      addGeneratorArgs(cmd, genNames);
      addInitialStateArgs(cmd);
      addParamArgs(cmd);

      Rivet::Commandline::addAnalysisArgs(cmd, analysesArg, analysesAllArg, anaNameConstraint);
      Rivet::Commandline::addHistoArgs(cmd, histoNameArg, histoTypeArg, histoTypeConstraint);
      Rivet::Commandline::addLoggingArgs(cmd, logsArg);

      TCLAP::ValueArg<size_t> 
        numEventsArg("n", "numevents", "Number of events to generate (10 by default)", 
                     false, 10, "num", cmd);

      TCLAP::ValueArg<string> 
        hepmcOutFileArg("o", "eventfile", "File to write HepMC events to (disabled by default)", 
                        false, "RivetGun.hepmc", "filename", cmd);
      
      TCLAP::SwitchArg disableRivetArg("R", "norivet", "Disable running of Rivet", cmd, false);

      // Parse the command line args and make a Configuration
      cmd.parse(argc, argv);

      // Populate the return value Configuration
      Configuration config;
      // Have to read params first in case run params are included
      useParamArgs(cmd, config);
      useGeneratorArgs(cmd, config);
      useInitialStateArgs(cmd, config);
      Rivet::Commandline::useAnalysisArgs(cmd, analysesArg, analysesAllArg, config.analyses);
      Rivet::Commandline::useHistoArgs(cmd, histoNameArg, histoTypeArg, config.histoName, config.histoFormat);
      Rivet::Commandline::useLoggingArgs(cmd, logsArg, config.logLevels);
      config.numEvents = numEventsArg.getValue();
      config.runRivet = ! disableRivetArg.getValue();
      if (hepmcOutFileArg.isSet()) {
        config.writeHepMC = true;
        config.hepmcOutFile = hepmcOutFileArg.getValue();
      }

      // Clean up the new'd arg parser objects
      discardArgObjects();

      // Add all the log levels from the command line into the logging framework
      Rivet::Log::setDefaultLevels(config.logLevels);

      return config;
    }


    void addGeneratorArgs(TCLAP::CmdLine& cmd, vector<string>& genNames) {
      genNameConstraint = new TCLAP::ValuesConstraint<string>(genNames);
      genNameArg = new TCLAP::ValueArg<string>("g", "generator", "Generator to use", true, "", genNameConstraint);
      genFileArg = new TCLAP::ValueArg<string>("G", "generatorfile", "HepML file defining the generator", true, "", "genfile");
      cmd.xorAdd(*genNameArg, *genFileArg);
      rngSeedArg = new TCLAP::ValueArg<size_t>("s", "seed", "random number generator seed", false, 271828, "seed", cmd);
    }


    void addInitialStateArgs(TCLAP::CmdLine& cmd) {
      vector<string> beamNames = Rivet::getKnownParticleNameNames();
      beamNameConstraint = new TCLAP::ValuesConstraint<string>(beamNames);
      beam1Arg = new TCLAP::ValueArg<string>("", "beam1", "Particle in beam 1 (PROTON by default)", false, "PROTON", beamNameConstraint, cmd);
      beam2Arg = new TCLAP::ValueArg<string>("", "beam2", "Particle in beam 2 (PROTON by default)", false, "PROTON", beamNameConstraint, cmd);
      mom1Arg = new TCLAP::ValueArg<double>("", "mom1", "Momentum of beam 1 in GeV (7000 by default)", false, 7000.0, "mom", cmd);
      mom2Arg = new TCLAP::ValueArg<double>("", "mom2", "Momentum of beam 2 in GeV (7000 by default)", false, 7000.0, "mom", cmd);
    }


    void addParamArgs(TCLAP::CmdLine& cmd) {
      paramsArg = new TCLAP::MultiArg<string>("p", "param", "Param in 'name:value' format", false, "param=value", cmd);
      paramsFilesArg = new TCLAP::MultiArg<string>("P", "paramsfile", "File containing parameter name, value pairs", false, "paramfile", cmd);
    }


    void useGeneratorArgs(TCLAP::CmdLine& cmd, Configuration& config) {
      if (!genNameArg || !genFileArg) return;
      if (genNameArg->isSet()) {
        config.generatorName = genNameArg->getValue();
      } else if (genFileArg->isSet()) {
        /// @todo Read HepML file
        throw runtime_error("HepML file reading is not yet supported. Sorry.");
      }
      if (rngSeedArg->isSet()) config.rngSeed = rngSeedArg->getValue();
    }


    void useInitialStateArgs(TCLAP::CmdLine& cmd, Configuration& config) {
      if (!beam1Arg || !mom1Arg || !beam2Arg || !mom2Arg) return;
      Rivet::ParticleNameMapR beamparticles = Rivet::getKnownParticleNamesR();
      try {
        if (beam1Arg->isSet()) config.beam1 = beamparticles[beam1Arg->getValue()];
        if (beam2Arg->isSet()) config.beam2 = beamparticles[beam2Arg->getValue()];
      } catch (exception& e) {
        throw runtime_error("Unknown beam particle: " + beam1Arg->getValue() + 
                            " or " + beam2Arg->getValue());
      }
      if (mom1Arg->isSet()) config.mom1 = mom1Arg->getValue();
      if (mom2Arg->isSet()) config.mom2 = mom2Arg->getValue();
    }

    
    void useParamArgs(TCLAP::CmdLine& cmd, Configuration& config) {
      if (!paramsArg || !paramsFilesArg) return;
      // Handle params from files/cin first...
      for (vector<string>::const_iterator pf = paramsFilesArg->getValue().begin(); 
           pf != paramsFilesArg->getValue().end(); ++pf) {
        if (*pf != "-") {
          /// @todo Logging at lower-than-default level
          // cout << "Using param file " << *pf << endl;
          std::ifstream in;
          in.open(pf->c_str());
          if (in.fail()) {
            throw runtime_error("Couldn't read from param file " + *pf);
          }
          handleConfigStream(in, config.params);
          in.close();
        } else {
          handleConfigStream(cin, config.params);
        }
      }
      // ...then overload the file params with specific command line params
      for (vector<string>::const_iterator p = paramsArg->getValue().begin(); 
           p != paramsArg->getValue().end(); ++p) {
        unsigned int breakpos = p->find("=");
        if (breakpos != string::npos) {
          string key = p->substr(0, breakpos);
          string value = p->substr(breakpos + 1, p->size() - breakpos - 1);
          config.params[key] = value;
        } else {
          throw runtime_error("Invalid parameter setting format: " + *p);
        }
      }

      // Handle parameters which should be intercepted by RivetGun
      for (map<string, string>::iterator p = config.params.begin();
           p != config.params.end(); ++p) {
        if (p->first.find("RG:") != string::npos) {
          //cout << "Found RG param: " << p->first << " = " << p->second << endl;
          stringstream pvalue(p->second);
          if (p->first == "RG:Mom1") {
            pvalue >> config.mom1;
          }
          if (p->first == "RG:Mom2") {
            pvalue >> config.mom2;
          }
          if (p->first == "RG:Beam1") {
            config.beam1 = Rivet::getKnownParticleNamesR()[p->second];
          }
          if (p->first == "RG:Beam2") {
            config.beam2 = Rivet::getKnownParticleNamesR()[p->second];
          }
          if (p->first == "RG:Seed") {
            pvalue >> config.rngSeed;
          }

          // Remove this param so it doesn't get passed to the generator
          config.params.erase(p);
        }
      }
    }


    void handleConfigStream(istream& in, map<string, string>& pmap) {
      while(in) {
        string line;
        getline(in, line);
        const string origline = line;
        //cout << "Input param line: " << line << endl;

        // Chop line after a # character
        size_t commentPosn = line.find_first_of("#");
        if (commentPosn != string::npos) {
          line = line.substr(0, commentPosn);
        }

        // Replace = signs with spaces
        while (line.find("=") != string::npos) {
          size_t eqpos = line.find("=");
          line.replace(eqpos, 1, " ");
        }

        // If effectively empty, then end
        if (line.find_first_not_of(' ') == string::npos) continue;

        // Get name and value and add to the param stack by finding the
        // punctuation tokens at either end of the string.
        const size_t first = line.find_first_not_of(' ');
        const size_t last = line.find_last_not_of(' ');
        line = line.substr(first, last-first+1);
        //cout << "Stripped line :" << line << ":" << endl;
        const size_t firstgap = line.find(' ');
        const size_t lastgap = line.rfind(' ');
        const string pname = line.substr(first, firstgap-first);
        const string pval = line.substr(lastgap+1, last-lastgap);
        //cout << ":" << pname << "=" << pval << ":" << endl;

        // Check if the remaining bit of string isn't emptiness
        const string remains = line.substr(firstgap, lastgap-firstgap+1);
        if (remains.find_first_not_of(' ') != string::npos) {
          cerr << remains << " - error at char #" << remains.find_first_not_of(' ') << endl;
          throw runtime_error("Param setting string '" + origline + "' is not in a valid format.");
        }

        // If all is okay, set the parameter
        if (!pname.empty() && !pval.empty()) {
          pmap[pname] = pval;
        } else {
          throw runtime_error("Param setting key or value is empty, somehow: '" + 
                              pname + "' = '" + pval + "'");
        }
      }
    }


  }  
}
