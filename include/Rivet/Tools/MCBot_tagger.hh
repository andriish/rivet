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
    bkg=0,
    top=1,
    V=2,
    H=3
  };


  //Class for doing MCBot tagging of reclustered Jets based on
  //https://gitlab.cern.ch/atlas-jetetmiss/tagging/mcbot
  //Class stores the model to ensure it is loaded from file
  //once and only once.
  //TODO: Is there a way to ensure model only reads once on multiple runs (e.g in GAMBIT)?
  // 
  class MCBot_tagger{
  public:

    //Constructor for string reference to json
    MCBot_tagger(const std::string& path_to_weights);

    //Constructor for string literal to json
    MCBot_tagger(const std::string&& path_to_weights);

    //Copy constructor
    //Will have same path to ref data, but initialises its own Neural Net
    MCBot_tagger(const MCBot_tagger& other);

    //Assignment operator - same principle as copy constructor
    MCBot_tagger operator=(const MCBot_tagger& tocopy);


    //Compute outputs for given inputs.
    //Its own function for historical reasons. Probably won need to be anymore
    void compute(const map<string, double>& inputs, map<string, double>& outputs) const;

    //Get scores from a Jet and its constituents
    void computeScores(const PseudoJet& totag, const Jets &constits, 
                                  std::map<string, double>& scoresOut) const;
    
    //Load a local copy of the network, compute output, and then throw away.
    //Useful for testing, probably useless long-term.
    void load_and_compute(const map<string, double>& inputs, map<string, double>& outputs) const;

    //Tag the provided PseudoJet
    //Should be const but I haven't gotten round to figuring out what to do with the 
    //threshold score map yet.
    MCBot_TagType tag(const PseudoJet& totag,const Jets &constituents);

    /// Get a logger object.
    Log& getLog() const;

  private:
    //the lwtnn neural net
    std::unique_ptr<lwt::LightweightNeuralNetwork> _lwg;

    //Path to the json that stores the network.
    std::string _path;
      
  };

  //Thresholds for individual scores.
  //TODO: Is there a "right" way to store these? Constexpr static or soemthing?
  // Score needs to be greater than to be tagged.
  std::map<string, double> _threshold = {{"PV", -0.2}, {"PH", 0.35}, {"Pt", 0.1}};

  //Thresholds for the double-tag tiebreaks.
  //A score greater than threshold indicates jet should be tagged as the second
  // of the pair of labels.
  std::map<string, double> _tiebreak_thresholds = {{"t_V", -0.3}, {"H_V", -0.55}, {"t_H", 0.2}};




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