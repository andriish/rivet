#include "Rivet/Rivet.hh"
#include "Rivet/HistoFormat.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/Configuration.hh"

#include "AGILe/Particle.hh"
#include "AGILe/Generator.hh"
#include "AGILe/Loader.hh"

#include <tclap/CmdLine.h>
using namespace TCLAP;

namespace {
  using namespace std;
  void handleConfigStream(istream& in, map<string, string>& pmap);
}


///////////////////////////////////////////


namespace Rivet {
  namespace Commandline {


    Configuration parseArgs(size_t argc, char** argv) {
      // The Configuration to be returned
      Configuration config;

      try {
        CmdLine cmd("Runs an event generator using the RivetGun and AGILe interface libraries", ' ', "1.0");

        // Add generator args
        vector<string> gens = AGILe::Loader::getAvailableGens();
        ValuesConstraint<string> genNameConstraint(gens);
        ValueArg<string> genNameArg("g", "generator", "Generator to use", true, "", &genNameConstraint);
        ValueArg<string> genFileArg("G", "generatorfile", "HepML file defining the generator", true, "", "genfile");
        cmd.xorAdd(genNameArg, genFileArg);
        ValueArg<size_t> rngSeedArg("s", "seed", "random number generator seed", false, 271828, "seed", cmd);
        ValueArg<string> eventFileArg("i", "ineventfile", "File containing HepMC events", false, "-", "filename", cmd);

        // Add initial state args
        vector<string> particles = Rivet::getParticleNames();
        ValuesConstraint<string> beamNameConstraint(particles);
        ValueArg<string> beam1Arg("", "beam1", "Particle in beam 1 (PROTON by default)", false, "PROTON", &beamNameConstraint, cmd);
        ValueArg<string> beam2Arg("", "beam2", "Particle in beam 2 (PROTON by default)", false, "PROTON", &beamNameConstraint, cmd);
        ValueArg<double> mom1Arg("", "mom1", "Momentum of beam 1 in GeV (7000 by default)", false, 7000.0, "mom", cmd);
        ValueArg<double> mom2Arg("", "mom2", "Momentum of beam 2 in GeV (7000 by default)", false, 7000.0, "mom", cmd);

        // Add generator param args
        MultiArg<string> paramsArg("p", "param", "Param in 'name:value' format", false, "param=value", cmd);
        MultiArg<string> paramsFilesArg("P", "paramsfile", "File containing parameter name, value pairs", false, "paramfile", cmd);

        // Add analysis args
        const set<string> tmp = Rivet::AnalysisLoader::getAllAnalysisNames();
        vector<string> knownAnalyses;
        for (set<string>::const_iterator i = tmp.begin(); i != tmp.end(); ++i) {
          knownAnalyses.push_back(*i);
          knownAnalyses.push_back("~" + *i);
        }
        ValuesConstraint<string> anaNameConstraint(knownAnalyses);
        MultiArg<string> analysesArg("a", "analysis", "A Rivet analysis to be run. Prefix name with a ~ to disable instead", false, &anaNameConstraint, cmd);
        SwitchArg analysesAllArg("A", "all_analyses", "Run all Rivet analyses (unless any are disabled)", cmd, false);

        // Add histogram args
        ValueArg<string> histoNameArg("H", "histoname", "Base name of histogram file, with no extension ('Rivet' by default)", false, "Rivet", "name", cmd);
        vector<string> histformats = Rivet::getKnownHistoFormatNames();
        ValuesConstraint<string> histoTypeConstraint(histformats);
        ValueArg<string> histoTypeArg("", "histotype", "Histogram output format (default is AIDA XML)", false, "AIDA", &histoTypeConstraint, cmd);

        // Add logging args
        const string mesg = "Set log level in 'name=level' format. The levels are INFO, DEBUG and WARNING";
        MultiArg<string> logsArg("l", "loglevel", mesg, false, "name=level", cmd);
        SwitchArg enableLogColorArg("", "color", "Use shell colors if possible (can be overridden by nocolor)", cmd, true);
        SwitchArg disableLogColorArg("", "nocolor", "Disable shell color escapes (useful for piping output to file)", cmd, false);

        // Add misc args
        ValueArg<size_t> numEventsArg("n", "numevents", "Number of events to generate (10 by default)", false, 10, "num", cmd);
        ValueArg<string> hepmcOutFileArg("o", "outeventfile", "File to write HepMC events to (disabled by default)", false, "RivetGun.hepmc", "filename", cmd);
        SwitchArg disableRivetArg("R", "norivet", "Disable running of Rivet", cmd, false);

        /////////////////////////////////////////////////////////


        // Parse the command line args and make a Configuration
        cmd.parse(argc, argv);


        /////////////////////////////////////////////////////////


        // Have to read params first in case run params are included
        // Handle params from files/cin first...
        for (vector<string>::const_iterator pf = paramsFilesArg.getValue().begin(); 
             pf != paramsFilesArg.getValue().end(); ++pf) {
          if (*pf != "-") {
            /// @todo Logging at lower-than-default level
            // cout << "Using param file " << *pf << endl;
            std::ifstream in;
            in.open(pf->c_str());
            if (in.fail() && pf->find("/") == string::npos) {
              // Try finding the file in the standard installed RivetGun share directory
              const string try2 = getRivetgunDataPath() + "/" + *pf;
              in.clear();
              in.open(try2.c_str());
              if (in.fail()) {
                throw runtime_error("Couldn't read from param file " + *pf + " or " + try2);
              }
            }
            handleConfigStream(in, config.params);
            in.close();
          } else {
            handleConfigStream(cin, config.params);
          }
        }
        // ...then overload the file params with specific command line params
        for (vector<string>::const_iterator p = paramsArg.getValue().begin(); 
             p != paramsArg.getValue().end(); ++p) {
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
              config.beam1 = Rivet::getParticleNameEnum(p->second);
            }
            if (p->first == "RG:Beam2") {
              config.beam2 = Rivet::getParticleNameEnum(p->second);
            }
            if (p->first == "RG:Seed") {
              pvalue >> config.rngSeed;
            }

            // Remove this param so it doesn't get passed to the generator
            config.params.erase(p);
          }
        }


        // Use generator args
        if (genNameArg.isSet()) {
          config.generatorName = genNameArg.getValue();
        } else if (genFileArg.isSet()) {
          /// @todo Read HepML file
          throw runtime_error("HepML file reading is not yet supported. Sorry.");
        }
        if (rngSeedArg.isSet()) config.rngSeed = rngSeedArg.getValue();
        if (eventFileArg.isSet()) {
          config.hepmcInFile = eventFileArg.getValue();
          config.readHepMC = true;
        }

        // Use initial state args
        try {
          if (beam1Arg.isSet()) config.beam1 = getParticleNameEnum(beam1Arg.getValue());
          if (beam2Arg.isSet()) config.beam2 = getParticleNameEnum(beam2Arg.getValue());
        } catch (exception& e) {
          throw runtime_error("Unknown beam particle: " + beam1Arg.getValue() + 
                              " or " + beam2Arg.getValue());
        }
        if (mom1Arg.isSet()) config.mom1 = mom1Arg.getValue();
        if (mom2Arg.isSet()) config.mom2 = mom2Arg.getValue();


        // Use analysis args
        // First handle the "enable all analyses" option...
        if (analysesAllArg.getValue()) {
          config.analyses = AnalysisLoader::getAllAnalysisNames();
        }
        // ...then handle individuals, including negations
        for (vector<string>::const_iterator ai = analysesArg.getValue().begin(); 
             ai != analysesArg.getValue().end(); ++ai) {
          string a = toUpper(*ai);
          try {
            // Check for analysis disabling with ~ prefix
            if (a.rfind("~", 0) == string::npos) {
              config.analyses.insert(a);
            } else {
              string aneg = a.substr(1, a.size()-1);
              config.analyses.erase(aneg);
            }
          } catch (std::exception& e) {
            throw std::runtime_error("Invalid analysis choice: " + *ai);
          }
        }


        // Use histo args
        // Get histogram filename
        config.histoName = histoNameArg.getValue();
        // Get histogram format
        if (histoTypeArg.getValue() == "AIDA") {
          config.histoFormat = AIDAML;
        } else if (histoTypeArg.getValue() == "FLAT") {
          config.histoFormat = FLAT;
        } else if (histoTypeArg.getValue() == "ROOT") {
          config.histoFormat = ROOT;
        }


        // Use logging args
        for (vector<string>::const_iterator l = logsArg.getValue().begin(); 
             l != logsArg.getValue().end(); ++l) {
          size_t breakpos = l->find("=");
          /// @todo Remove backwards compatible ":" delimiter
          if (breakpos == string::npos) breakpos = l->find(":");
          if (breakpos != string::npos) {
            string key = l->substr(0, breakpos);
            string value = l->substr(breakpos + 1, l->size() - breakpos - 1);
            config.logLevels[key] = Log::getLevelFromName(value);
          } else {
            throw runtime_error("Invalid log setting format: " + *l);
          }
        }
        if (enableLogColorArg.getValue()) config.useLogColors = true;
        if (disableLogColorArg.getValue()) config.useLogColors = false;

        config.numEvents = numEventsArg.getValue();
        config.runRivet = ! disableRivetArg.getValue();
        if (hepmcOutFileArg.isSet()) {
          config.writeHepMC = true;
          config.hepmcOutFile = hepmcOutFileArg.getValue();
        }
      } catch (const TCLAP::ArgException& e) { 
        stringstream ss;
        ss << "Command line error: " << e.error() << " for arg " << e.argId(); 
        throw runtime_error(ss.str());
      }

      return config;
    }


  }  
}


//////////////////////////////////////////////


namespace {
  using namespace std;
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
