#include "Rivet/Tools/MCBot_tagger.hh"
#include "lwtnn/NanReplacer.hh"

namespace Rivet{

  

  //Constructor for string reference to json
  MCBot_tagger::MCBot_tagger(const std::string& path_to_weights) : _path(path_to_weights){
    std::ifstream input(_path);
    auto config = lwt::parse_json(input);
    input.close();

    _lwg = std::make_unique<lwt::LightweightNeuralNetwork>(config.inputs, config.layers, config.outputs);
    lwt::NanReplacer replacer(config.defaults, lwt::rep::all);
  }

  //Constructor for string literal to json
  MCBot_tagger::MCBot_tagger(const std::string&& path_to_weights) : _path(std::move(path_to_weights)){
    std::ifstream input(_path);
    auto config = lwt::parse_json(input);
    input.close();

    _lwg = std::make_unique<lwt::LightweightNeuralNetwork>(config.inputs, config.layers, config.outputs);
    lwt::NanReplacer replacer(config.defaults, lwt::rep::all);
  }

  //Copy constructor
  //Will have same path to ref data, but initialises its own Neural Net
  MCBot_tagger::MCBot_tagger(const MCBot_tagger& other) : MCBot_tagger(other._path){
  }

  //Assignment operator - same principle as copy constructor
  MCBot_tagger MCBot_tagger::operator=(const MCBot_tagger& other){
    return MCBot_tagger(other._path);
  }


  void MCBot_tagger::load_and_compute(map<string, double>& inputs, map<string, double>& outputs) const{

    //Is it insanely dirty to load every time? yes.
    //But I need to get the blasted thing working and we can iron out details later.
    std::ifstream input(_path);
    auto config = lwt::parse_json(input);
    input.close();
    lwt::LightweightNeuralNetwork lwg(config.inputs, config.layers, config.outputs);
    lwt::NanReplacer replacer(config.defaults, lwt::rep::all);

    outputs = lwg.compute(inputs);

  }


  void MCBot_tagger::compute(map<string, double>& inputs, map<string, double>& outputs) const{
    outputs = _lwg->compute(inputs);
  }


  MCBot_TagType MCBot_tagger::tag(const PseudoJet& totag, const Jets &constits){
    Jets leadingsubjets = sortByPt(constits);

    if (leadingsubjets.size() < 3){
      MSG_WARNING("Unable to tag vRC jet with only " << leadingsubjets.size() << " consituents!");
      return MCBot_TagType::bkg;
    }

    std::map<string, double> input_vals = {
      {"rcjet_pt", totag.pt()},
      {"rcjet_numConstituents", static_cast<double>(constits.size())},
      {"rcjet_m", totag.m()},
      
      {"sjet_1_mv2c10_binned", static_cast<double>(hasBTag()(leadingsubjets[2]))},//todo: double check this
      {"sjet_1_e", leadingsubjets[0].E()},
      {"sjet_1_phi", leadingsubjets[0].phi()},
      {"sjet_1_eta", leadingsubjets[0].eta()},
      {"sjet_1_pt", leadingsubjets[0].pt()},

      {"sjet_2_mv2c10_binned", static_cast<double>(hasBTag()(leadingsubjets[2]))},//todo: double check this
      {"sjet_2_e", leadingsubjets[1].E()},
      {"sjet_2_phi", leadingsubjets[1].phi()},
      {"sjet_2_eta", leadingsubjets[1].eta()},
      {"sjet_2_pt", leadingsubjets[1].pt()},

      {"sjet_3_mv2c10_binned", static_cast<double>(hasBTag()(leadingsubjets[2]))},//todo: double check this
      {"sjet_3_e", leadingsubjets[2].E()},
      {"sjet_3_phi", leadingsubjets[2].phi()},
      {"sjet_3_eta", leadingsubjets[2].eta()},
      {"sjet_3_pt", leadingsubjets[2].pt()}
    };


    std::map<string, double> outputs;
    //load_and_compute(input_vals, outputs);
    compute(input_vals, outputs);


    double PV=log10(outputs["dnnOutput_V"]/
      (0.9*outputs["dnnOutput_light"]+0.05*outputs["dnnOutput_top"]+0.05*outputs["dnnOutput_H"]));
    double PH=log10(outputs["dnnOutput_H"]/
      (0.9*outputs["dnnOutput_light"]+0.05*outputs["dnnOutput_top"]+0.05*outputs["dnnOutput_V"]));
    double Ptop=log10(outputs["dnnOutput_top"]/
      (0.9*outputs["dnnOutput_light"]+0.05*outputs["dnnOutput_V"]+0.05*outputs["dnnOutput_H"]));

    //Are values above the threshold value?
    bool isV = (PV > _threshold["PV"]);
    bool isH = (PH > _threshold["PH"]);
    bool istop = (Ptop > _threshold["Ptop"]);

    MSG_INFO("(DV, DH, Dtop, Dlight) = (" << outputs["dnnOutput_V"] << ", " << 
        outputs["dnnOutput_H"] << ", " << outputs["dnnOutput_top"] << ", " <<
         outputs["dnnOutput_light"]  << ")");
    MSG_INFO("(PV, PH, Ptop) = (" << PV << ", " << PH << ", " << Ptop << ")");
    MSG_INFO("(isV, isH, istop) = (" << isV << ", " << isH << ", " << istop << ")");

    //Return values
    if ((!isV) && (!isH) && (!istop)){
      return MCBot_TagType::bkg;
    }
    else if (isV && (!isH) && (!istop)){
      return MCBot_TagType::V;
    }
    //Note triple tagged is counted as Higgs
    else if (((!isV) && isH && (!istop)) || (isV && isH && istop)){
      return MCBot_TagType::H;
    }
    else if ((!isH) && (!isV) && istop){
      return MCBot_TagType::top;
    }
    //Double tag cases require a tiebreak.
    else if (isH && isV){
      double tiebreakScore = log10(outputs["dnnOutput_V"]/outputs["dnnOutput_H"]);
      MSG_INFO("Tie break score: " << tiebreakScore << " (threshold " << _tiebreak_thresholds["H_V"] << ")");

      return (tiebreakScore > _tiebreak_thresholds["H_V"]) ? MCBot_TagType::V : MCBot_TagType::H;
    }
    else if (isV && istop){
      double tiebreakScore = log10(outputs["dnnOutput_V"]/outputs["dnnOutput_top"]);
      return (tiebreakScore > _tiebreak_thresholds["t_V"]) ? MCBot_TagType::V : MCBot_TagType::top;
    }
    else {
      double tiebreakScore = log10(outputs["dnnOutput_H"]/outputs["dnnOutput_top"]);
      return (tiebreakScore > _tiebreak_thresholds["t_H"]) ? MCBot_TagType::H : MCBot_TagType::top;
    }
  }

  Log& MCBot_tagger::getLog() const {
    return Rivet::Log::getLog("Rivet.MCBot_tagger");
  }





}    