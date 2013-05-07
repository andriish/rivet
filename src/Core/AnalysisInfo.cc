#include "Rivet/Rivet.hh"
#include "Rivet/RivetBoost.hh"
#include "Rivet/AnalysisInfo.hh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/Tools/RivetPaths.hh"
#include "Rivet/Tools/Logging.hh"
#include "yaml-cpp/yaml.h"
#include <iostream>
#include <fstream>
#include <unistd.h>

namespace Rivet {


  namespace {
    Log& getLog() {
      return Log::getLog("Rivet.AnalysisInfo");
    }
  }


  /// Static factory method
  AnalysisInfo* AnalysisInfo::make(const std::string& ananame) {
    // Returned AI, in semi-null state
    AnalysisInfo* ai = new AnalysisInfo();
    ai->_beams += make_pair(PID::ANY, PID::ANY);
    ai->_name = ananame;

    /// If no ana data file found, return null AI
    const string datapath = findAnalysisInfoFile(ananame + ".info");
    if (datapath.empty()) {
      MSG_DEBUG("No datafile " << ananame + ".info found");
      return ai;
    }

    // Read data from YAML document
    MSG_DEBUG("Reading analysis data from " << datapath);
    YAML::Node doc;
    try {
      doc = YAML::LoadFile(datapath);
    } catch (const YAML::ParserException& ex) {
      MSG_ERROR("Parse error when reading analysis data from " << datapath << " (" << ex.what() << ")");
      return ai;
    }

    #define THROW_INFOERR(KEY) throw InfoError("Problem in info parsing while accessing key " + string(KEY) + " in file " + datapath)
    #define TRY_GETINFO(KEY, VAR) try { if (doc[KEY] && !doc[KEY].IsNull()) ai->_ ## VAR = doc[KEY].as<string>(); } catch (...) { THROW_INFOERR(KEY); }

    TRY_GETINFO("Name", name);
    // if (doc["Summary"]) ai->_summary = doc["Summary"].as<string>();
    TRY_GETINFO("Summary", summary);
    // if (doc["Status"]) ai->_status = doc["Status"].as<string>();
    TRY_GETINFO("Status", status);
    // if (doc["RunInfo"]) ai->_runInfo = doc["RunInfo"].as<string>();
    TRY_GETINFO("RunInfo", runInfo);
    // if (doc["Description"]) ai->_description = doc["Description"].as<string>();
    TRY_GETINFO("Description", description);
    // if (doc["Experiment"]) ai->_experiment = doc["Experiment"].as<string>();
    TRY_GETINFO("Experiment", experiment);
    // if (doc["Collider"]) ai->_collider = doc["Collider"].as<string>();
    TRY_GETINFO("Collider", collider);
    // if (doc["Year"]) ai->_year = doc["Year"].as<string>();
    TRY_GETINFO("Year", year);
    // if (doc["SpiresID"]) ai->_spiresId = doc["SpiresID"].as<string>();
    TRY_GETINFO("SpiresID", spiresId);
    // if (doc["InspireID"]) ai->_inspireId = doc["InspireID"].as<string>();
    TRY_GETINFO("InspireID", inspireId);
    // if (doc["BibKey"]) ai->_bibKey = doc["BibKey"].as<string>();
    TRY_GETINFO("BibKey", bibKey);
    // if (doc["BibTeX"]) ai->_bibTeX = doc["BibTeX"].as<string>();
    TRY_GETINFO("BibTeX", bibTeX);

    try {
      if (doc["NeedsCrossSection"]) ai->_needsCrossSection = doc["NeedsCrossSection"].as<bool>();
      else if (doc["NeedsCrossSection"]) ai->_needsCrossSection = doc["NeedCrossSection"].as<bool>();
    } catch (...) {
      THROW_INFOERR("NeedsCrossSection|NeedCrossSection");
    }

    try {
      if (doc["Authors"]) {
        const YAML::Node& authors = doc["Authors"];
        for (size_t i = 0; i < authors.size(); ++i) ai->_authors += authors[i].as<string>();
      }
    } catch (...) { THROW_INFOERR("Authors"); }

    try {
      if (doc["References"]) {
        const YAML::Node& refs = doc["References"];
        for (size_t i = 0; i < refs.size(); ++i) ai->_references += refs[i].as<string>();
      }
    } catch (...) { THROW_INFOERR("References"); }

    try {
      if (doc["ToDo"]) {
        const YAML::Node& todos = doc["ToDo"];
        for (size_t i = 0; i < todos.size(); ++i) ai->_todos += todos[i].as<string>();
      }
    } catch (...) { THROW_INFOERR("ToDo"); }

    try {
      if (doc["Beams"]) {
        const YAML::Node& beams = doc["Beams"];
        vector<PdgIdPair> beam_pairs;
        if (beams.size() == 2 && beams[0].IsScalar() && beams[0].IsScalar()) {
          beam_pairs += PID::make_pdgid_pair(beams[0].as<string>(), beams[1].as<string>());
        } else {
          for (size_t i = 0; i < beams.size(); ++i) {
            const YAML::Node& bp = beams[i];
            if (bp.size() != 2 || !bp[0].IsScalar() || !bp[0].IsScalar())
              throw InfoError("Beam ID pairs have to be either a 2-tuple or a list of 2-tuples of particle names");
            beam_pairs += PID::make_pdgid_pair(bp[0].as<string>(), bp[1].as<string>());
          }
        }
        ai->_beams = beam_pairs;
      }
    } catch (...) { THROW_INFOERR("beams"); }

    try {
      if (doc["Energies"]) {
        vector< pair<double,double> > beam_energy_pairs;
        for (size_t i = 0; i < doc["Energies"].size(); ++i) {
          const YAML::Node& be = doc["Energies"][i];
          if (be.IsScalar()) {
            // If beam energy is a scalar, then assume symmetric beams each with half that energy
            beam_energy_pairs += make_pair(be.as<double>()/2.0, be.as<double>()/2.0);
          } else if (be.IsSequence()) {
            if (be.size() != 2)
              throw InfoError("Beam energies have to be a list of either numbers or pairs of numbers");
            beam_energy_pairs += make_pair(be[0].as<double>(), be[1].as<double>());
          } else {
            throw InfoError("Beam energies have to be a list of either numbers or pairs of numbers");
          }
        }
        ai->_energies = beam_energy_pairs;
      }
    } catch (...) { THROW_INFOERR("Energies"); }

    #undef TRY_GETINFO
    #undef THROW_INFOERR

    MSG_TRACE("AnalysisInfo pointer = " << ai);
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
