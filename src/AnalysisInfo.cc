#include "Rivet/Rivet.hh"
#include "Rivet/AnalysisInfo.hh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/Tools/Logging.hh"
#include "yaml.h"
#include <iostream>
#include <fstream>
#include <unistd.h>

namespace Rivet {


  /// Ideas: 
  ///  * search RIVET_DATA_PATH etc. for <name>.info.yaml
  ///  * how to determine the name?
  ///  * only populate pointer on Analysis when requested
  ///  * use smart pointer: deletes automatically when Analysis 
  ///    goes out of scope


  /// Static factory method
  AnalysisInfo* AnalysisInfo::make(const std::string& name) {
    // Search metadata path and read first matching file
    string datapath = getRivetDataPath() + "/" + name + ".info";
    Log::getLog("Rivet.Analysis") 
      << Log::DEBUG << "Reading analysis data from " << datapath << endl;
    if (access(datapath.c_str(), R_OK) != 0) return 0;

    // Read data from YAML document
    std::ifstream io(datapath.c_str());
    YAML::Parser parser(io);
    YAML::Node doc;
    try {
      parser.GetNextDocument(doc);
      //cout << doc << endl;
    } catch (const YAML::ParserException& ex) {
      Log::getLog("Rivet.Analysis") 
        << Log::ERROR << "Parse error when reading analysis data from " 
        << datapath << endl;
      return 0;
    }

    AnalysisInfo* ai = new AnalysisInfo();
    for (YAML::Iterator it = doc.begin(); it != doc.end(); ++it) {
      string key;
      it.first() >> key;
      Log::getLog("Rivet.Analysis") 
        << Log::TRACE << key << ": " << it.second() << endl;
      try {
        if (key == "Name") {
          it.second() >> ai->_name;
        } else if (key == "Summary") {
          it.second() >> ai->_summary;
        } else if (key == "Experiment") {
          it.second() >> ai->_experiment;
        } else if (key == "Beams") {
          // it.second() >> ai->_beams;
        } else if (key == "NeedCrossSection") {
          // it.second() >> ai->_needsCrossSection;
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
    //cout << *ai << endl;
    return ai;
  }


  string toString(const AnalysisInfo& ai) {
    stringstream ss;
    ss << ai.name();
    ss << " - " << ai.summary();
    ss << " (" << ai.status() << ")";
    return ss.str();
  }


}
