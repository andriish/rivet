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
      #if YAMLCPP_API_VERSION == 3
      std::ifstream file(datapath);
      doc = YAML::Load(file);
      #elif YAMLCPP_API_VERSION == 5
      doc = YAML::LoadFile(datapath);
      #endif
    } catch (const YAML::ParserException& ex) {
      MSG_ERROR("Parse error when reading analysis data from " << datapath << " (" << ex.what() << ")");
      return ai;
    }

    #define THROW_INFOERR(KEY) throw InfoError("Problem in info parsing while accessing key " + string(KEY) + " in file " + datapath)

    // Simple scalars (test for nullness before casting)
    #if YAMLCPP_API_VERSION == 3
    /// @todo Fix
    #define TRY_GETINFO(KEY, VAR) try { if (doc.FindValue(KEY) { string val; doc[KEY] >> val; ai->_ ## VAR = val; } catch (...) { THROW_INFOERR(KEY); }
    #elif YAMLCPP_API_VERSION == 5
    #define TRY_GETINFO(KEY, VAR) try { if (doc[KEY] && !doc[KEY].IsNull()) ai->_ ## VAR = doc[KEY].as<string>(); } catch (...) { THROW_INFOERR(KEY); }
    #endif
    TRY_GETINFO("Name", name);
    TRY_GETINFO("Summary", summary);
    TRY_GETINFO("Status", status);
    TRY_GETINFO("RunInfo", runInfo);
    TRY_GETINFO("Description", description);
    TRY_GETINFO("Experiment", experiment);
    TRY_GETINFO("Collider", collider);
    TRY_GETINFO("Year", year);
    TRY_GETINFO("SpiresID", spiresId);
    TRY_GETINFO("InspireID", inspireId);
    TRY_GETINFO("BibKey", bibKey);
    TRY_GETINFO("BibTeX", bibTeX);
    #undef TRY_GETINFO

    // Sequences (test the seq *and* each entry for nullness before casting)
    #if YAMLCPP_API_VERSION == 3
    /// @todo Fix
    #define TRY_GETINFO_SEQ(KEY, VAR) try { \
        if (const YAML::Node* VAR = doc.FindValue(KEY)) {                                       \
          for (size_t i = 0; i < VAR.size(); ++i) {                     \
            string val; VAR[i] >> val; ai->_ ## VAR += val;             \
          } } } catch (...) { THROW_INFOERR(KEY); }
    #elif YAMLCPP_API_VERSION == 5
    #define TRY_GETINFO_SEQ(KEY, VAR) try { \
        if (doc[KEY] && !doc[KEY].IsNull()) {                           \
          const YAML::Node& VAR = doc[KEY];                             \
          for (size_t i = 0; i < VAR.size(); ++i)                       \
            if (!VAR[i].IsNull()) ai->_ ## VAR += VAR[i].as<string>();  \
        } } catch (...) { THROW_INFOERR(KEY); }
    #endif
    TRY_GETINFO_SEQ("Authors", authors);
    TRY_GETINFO_SEQ("References", references);
    TRY_GETINFO_SEQ("ToDo", todos);
    #undef TRY_GETINFO_SEQ


    // A boolean with some name flexibility
    try {
      #if YAMLCPP_API_VERSION == 3
      /// @todo Fix
      #elif YAMLCPP_API_VERSION == 5
      if (doc["NeedsCrossSection"]) ai->_needsCrossSection = doc["NeedsCrossSection"].as<bool>();
      else if (doc["NeedsCrossSection"]) ai->_needsCrossSection = doc["NeedCrossSection"].as<bool>();
      #endif
    } catch (...) {
      THROW_INFOERR("NeedsCrossSection|NeedCrossSection");
    }


    // Beam particle identities
    try {
      #if YAMLCPP_API_VERSION == 3

      /// @todo Fix

      #elif YAMLCPP_API_VERSION == 5

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

      #endif
    } catch (...) { THROW_INFOERR("Beams"); }


    // Beam energies
    try {
      #if YAMLCPP_API_VERSION == 3

      /// @todo Fix

      #elif YAMLCPP_API_VERSION == 5

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

      #endif

    } catch (...) { THROW_INFOERR("Energies"); }

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
