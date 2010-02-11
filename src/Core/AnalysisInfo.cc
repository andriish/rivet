#include "Rivet/Rivet.hh"
#include "Rivet/AnalysisInfo.hh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/Tools/Logging.hh"
#include "yaml-cpp/yaml.h"
#include <iostream>
#include <fstream>
#include <unistd.h>

namespace Rivet {


  /// Ideas:
  ///  * search RIVET_INFO_PATH etc. for <name>.info.yaml
  ///  * how to determine the name?
  ///  * only populate pointer on Analysis when requested
  ///  * use smart pointer: deletes automatically when Analysis
  ///    goes out of scope


  /// Static factory method
  AnalysisInfo* AnalysisInfo::make(const std::string& ananame) {
    // Build the list of directories to search
    vector<string> dirs;
    char* env = 0;
    // First try to use the Rivet data path variable
    env = getenv("RIVET_INFO_PATH");
    if (env) dirs += split(env);
    // Then try to use the Rivet data install path
    dirs += getRivetDataPath();
    // And the current dir
    dirs += ".";

    bool found = false;
    string datapath = "";
    foreach (const string& d, dirs) {
      if (d.empty()) continue;
      /// @todo Use system-independent separator (e.g. Boost.File)
      datapath = d + "/" + ananame + ".info";
      Log::getLog("Rivet.AnalysisInfo")
        << Log::TRACE << "Looking for analysis data file '" << datapath << "'" << endl;
      if (access(datapath.c_str(), R_OK) == 0) {
        found = true;
        break;
      }
    }

    // Returned AI, in semi-null state
    AnalysisInfo* ai = new AnalysisInfo();
    ai->_beams = make_pair(ANY,ANY);
    ai->_name = ananame;

    /// If no ana data file found, return null AI
    if (!found) return ai;

    // Read data from YAML document
    Log::getLog("Rivet.AnalysisInfo")
      << Log::DEBUG << "Reading analysis data from " << datapath << endl;
    std::ifstream io(datapath.c_str());
    YAML::Parser parser(io);
    YAML::Node doc;
    try {
      parser.GetNextDocument(doc);
      //cout << doc << endl;
    } catch (const YAML::ParserException& ex) {
      Log::getLog("Rivet.AnalysisInfo")
        << Log::ERROR << "Parse error when reading analysis data from "
        << datapath << endl;
      return ai;
    }

    for (YAML::Iterator it = doc.begin(); it != doc.end(); ++it) {
      string key;
      it.first() >> key;
      stringstream sec;
      sec << it.second();
      const string secstr = sec.str().substr(0, sec.str().length()-1);
      Log::getLog("Rivet.AnalysisInfo")
        << Log::TRACE << key << ": " << secstr << endl;
      try {
        if (key == "Name") {
          it.second() >> ai->_name;
        } else if (key == "Summary") {
          it.second() >> ai->_summary;
        } else if (key == "Experiment") {
          it.second() >> ai->_experiment;
        } else if (key == "Beams") {
          const YAML::Node& beams = it.second();
          vector<ParticleName> beampair;
          for (YAML::Iterator b = beams.begin(); b != beams.end(); ++b) {
            string bstr;
            *b >> bstr;
            ParticleName beamname = getParticleNameEnum(bstr);
            beampair += beamname;
          }
          assert(beampair.size() == 2);
          ai->_beams = make_pair<ParticleName,ParticleName>(beampair[0], beampair[1]);
        // } else if (key == "NeedCrossSection") {
        //   // it.second() >> ai->_needsCrossSection;
        } else if (key == "Energies") {
          const YAML::Node& energies = it.second();
          vector<pair<double,double> > beam_energy_pairs;
          for (YAML::Iterator be = energies.begin(); be != energies.end(); ++be) {
            if (be->GetType() == YAML::CT_SCALAR) {
              // If beam energy is a scalar, then assume symmetric beams each with half that energy
              double sqrts;
              *be >> sqrts;
              beam_energy_pairs += make_pair(sqrts/2.0, sqrts/2.0);
            } else if (be->GetType() == YAML::CT_SEQUENCE) {
              const YAML::Node& beamenergy_strs = be.second();
              vector<double> beamenergies;
              for (YAML::Iterator e = beamenergy_strs.begin(); e != beamenergy_strs.end(); ++e) {
                double beamenergy;
                *e >> beamenergy;
                beamenergies += beamenergy;
              }
              assert(beamenergies.size() == 2);
              beam_energy_pairs += make_pair(beamenergies[0], beamenergies[1]);
            } else {
              assert(0 && "Beam energies have to be a list of either numbers or pairs of numbers");
            }
          }
          ai->_energies = beam_energy_pairs;
        } else if (key == "Collider") {
          it.second() >> ai->_collider;
        } else if (key == "SpiresID") {
          it.second() >> ai->_spiresId;
        } else if (key == "Status") {
          it.second() >> ai->_status;
        } else if (key == "RunInfo") {
          it.second() >> ai->_runInfo;
        } else if (key == "Description") {
          it.second() >> ai->_description;
        } else if (key == "Year") {
          it.second() >> ai->_year;
        } else if (key == "Authors") {
          const YAML::Node& authors = it.second();
          for (YAML::Iterator a = authors.begin(); a != authors.end(); ++a) {
            string astr;
            *a >> astr;
            ai->_authors += astr;
          }
        } else if (key == "References") {
          const YAML::Node& refs = it.second();
          for (YAML::Iterator r = refs.begin(); r != refs.end(); ++r) {
            string rstr;
            *r >> rstr;
            ai->_references += rstr;
          }
        }
      } catch (const YAML::RepresentationException& ex) {
        Log::getLog("Rivet.Analysis")
          << Log::WARN << "Type error when reading analysis data '"
          << key << "' from " << datapath << endl;
      }
    }
    Log::getLog("Rivet.AnalysisInfo") << Log::DEBUG << ai << endl;
    return ai;
  }


  string toString(const AnalysisInfo& ai) {
    stringstream ss;
    ss << ai.name();
    ss << " - " << ai.summary();
    // ss << " - " << ai.beams();
    // ss << " - " << ai.energies();
    ss << " (" << ai.status() << ")";
    return ss.str();
  }


}
