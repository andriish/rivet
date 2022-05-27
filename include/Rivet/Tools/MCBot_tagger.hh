#ifndef RIVET_MCBOT_TAGGER_HH
#define RIVET_MCBOT_TAGGER_HH

//TODO: Eventually perhpas belongs in jetnutils?

#include <fstream>

#include "Rivet/Jet.hh"
#include "Rivet/Tools/JetUtils.hh"


#include "lwtnn/LightweightNeuralNetwork.hh"
#include "lwtnn/parse_json.hh"

namespace Rivet {

  //enum class for the four possible output classification types of the 
  // MCBot tagger
  enum class MCBot_TagType{
    bkg,
    top,
    V,
    H
  };


  //Class for doing MCBot tagging of reclustered Jets based on
  //https://gitlab.cern.ch/atlas-jetetmiss/tagging/mcbot
  //Class stores the model to ensure it is loaded from file
  //once and only once.
  //TODO: Is there a way to ensure model only reads once on multiple runs (e.g in GAMBIT)?
  // 
  class MCBot_tagger{
  public:
    //Constructor
    MCBot_tagger(const std::string& path_to_weights);

    MCBot_tagger(const std::string&& path_to_weights);

    void load_and_compute(map<string, double>& inputs, map<string, double>& outputs);

    //Tag the provided PseudoJet
    MCBot_TagType tag(const PseudoJet& totag,const Jets &constituents);

  private:
    //the lwtnn neural net
    //std::shared_ptr<lwt::LightweightNeuralNetwork> _lwg;

    bool _lwNNloaded;
    std::string _path;

    //Thresholds for individual scores.
    // Score needs to be greater than to be tagged.
    std::map<string, double> _threshold = {{"PV", -0.2}, {"PH", 0.35}, {"Pt", 0.1}};

    //Thresholds for the double-tag tiebreaks.
    //A score greater than threshold indicates jet should be tagged as the second
    // of the pair of labels.
    std::map<string, double> _tiebreak_thresholds ={{"t_V", -0.3}, {"H_V", -0.55}, {"t_H", 0.2}};


    /// Get a logger object.
    Log& getLog() const;
      
  };


  // struct HasMCBot_VTag : BoolJetFunctor {
  //   HasMCBot_VTag(MCBot_tagger& tagger){
  //     _tagger = tagger;
  //   }

  //   bool operator(const)

  // private:
  //   MCBot_tagger& _tagger

  // };
}

#endif