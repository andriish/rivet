#include "Rivet/Tools/MCBot_tagger.hh"

namespace Rivet{

  MCBot_tagger::MCBot_tagger(std::string& path_to_weights){
    std::ifstream input(path_to_weights);
    auto config = lwt::parse_json(input);
    _lwg = std::make_shared<lwt::LightweightNeuralNetwork>(config.inputs, config.layers, config.outputs);
  }

  MCBot_tagger::MCBot_tagger(std::string path_to_weights){
    std::cout << __FILE__ << ": " << __LINE__ << std::endl;
    std::ifstream input(path_to_weights);
    std::cout << __FILE__ << ": " << __LINE__ << std::endl;
    auto config = lwt::parse_json(input);
    std::cout << __FILE__ << ": " << __LINE__ << std::endl;
    _lwg = std::make_shared<lwt::LightweightNeuralNetwork>(config.inputs, config.layers, config.outputs);
    std::cout << __FILE__ << ": " << __LINE__ << std::endl;
  }


  MCBot_TagType MCBot_tagger::tag(const PseudoJet& totag,const Jets &constituents){
    std::cout << __FILE__ << ": " << __LINE__ << std::endl;
    Jets leadingsubjets = sortByPt(constituents);

    std::map<string, double> inputs = {
      {"rcjet_pt", totag.pt()},
      {"rcjet_numConstituents", static_cast<double>(totag.constituents().size())},
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
      {"sjet_4_pt", leadingsubjets[2].pt()}
    };

    std::map<string, double> outputs = _lwg->compute(inputs);

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

    //Return values
    if (!isV && !isH && !istop){
      return MCBot_TagType::bkg;
    }
    else if (isV && !isH && !istop){
      return MCBot_TagType::V;
    }
    //Note triple tagged is counted as Higgs
    else if ((!isV && isH && !istop) || (isV && isH && istop)){
      return MCBot_TagType::H;
    }
    else if (!isH && !isV && istop){
      return MCBot_TagType::top;
    }
    //Double tag cases require a tiebreak.
    else if (isH && isV){
      double tiebreakScore = log10(outputs["dnnOutput_V"]/outputs["dnnOutput_H"]);
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

    std::cout << __FILE__ << ": " << __LINE__ << std::endl;
  }


}    