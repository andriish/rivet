#include "Rivet/Tools/MCBot_tagger.hh"
#include "lwtnn/NanReplacer.hh"

namespace Rivet{

  

  
  MCBot_tagger::MCBot_tagger(const std::string& path_to_weights) : _path(path_to_weights){
    //Use a try catch block because this can fail (e.g. bad path, wring LD_LIBRARY_PATH for lwtnn)
    // and if it does fail, the error messages are always very misleading.
    try {
      std::ifstream input(_path);
      auto config = lwt::parse_json(input);
      input.close();

      _lwg = std::make_unique<lwt::LightweightNeuralNetwork>(config.inputs, config.layers, config.outputs);
      //TODO: I'm not sure what the point of the NaN replacer is, but I'm using it to copy 
      //the original analyses work flow. 
      lwt::NanReplacer replacer(config.defaults, lwt::rep::all);
    }
    catch (Exception &e){
      throw "Error in initialising MCBot tagger. Does the json file exist, and is it correctly structured?";
    }
  }

  MCBot_tagger::MCBot_tagger(const std::string&& path_to_weights) : _path(std::move(path_to_weights)){
  //Use a try catch block because this can fail (e.g. bad path, wring LD_LIBRARY_PATH for lwtnn)
    // and if it does fail, the error messages are always very misleading.
    try {
    std::ifstream input(_path);
    auto config = lwt::parse_json(input);
    input.close();

    _lwg = std::make_unique<lwt::LightweightNeuralNetwork>(config.inputs, config.layers, config.outputs);
    //TODO: I'm not sure what the point of the NaN replacer is, but I'm using it to copy 
      //the original analyses work flow. 
    lwt::NanReplacer replacer(config.defaults, lwt::rep::all);
    }
    catch (Exception &e){
      throw "Error in initialising MCBot tagger. Does the json file exist, and is it correctly structured?";
    }
  }

  MCBot_tagger::MCBot_tagger(const MCBot_tagger& other) : MCBot_tagger(other._path){
  }

  MCBot_tagger MCBot_tagger::operator=(const MCBot_tagger& other){
    return MCBot_tagger(other._path);
  }


  void MCBot_tagger::load_and_compute(const map<string, double>& inputs, map<string, double>& outputs) const{

    //Is it insanely dirty to load every time? yes.
    //But I need to get the blasted thing working and we can iron out details later.
    std::ifstream input(_path);
    auto config = lwt::parse_json(input);
    input.close();
    lwt::LightweightNeuralNetwork lwg(config.inputs, config.layers, config.outputs);
    lwt::NanReplacer replacer(config.defaults, lwt::rep::all);

    outputs = lwg.compute(inputs);

  }


  void MCBot_tagger::compute(const map<string, double>& inputs, map<string, double>& outputs) const{
    outputs = _lwg->compute(inputs);
  }


  void MCBot_tagger::computeScores(const PseudoJet& totag, const Jets &constits, 
                                  std::map<string, double>& scoresOut) const{

    Jets leadingsubjets = sortByPt(constits);

    //MCBot is trained in MeV units (I think)
    std::map<string, double> input_vals = {
      {"rcjet_pt", totag.pt()/MeV},
      {"rcjet_numConstituents", static_cast<double>(constits.size())},
      {"rcjet_m", totag.m()/MeV},
      
      {"sjet_1_mv2c10_binned", static_cast<double>(hasBTag()(leadingsubjets[0]))},
      {"sjet_1_e", leadingsubjets[0].E()/MeV},
      {"sjet_1_phi", leadingsubjets[0].phi()},
      {"sjet_1_eta", leadingsubjets[0].eta()},
      {"sjet_1_pt", leadingsubjets[0].pt()/MeV}
    };

    if (leadingsubjets.size() >= 2){
      input_vals["sjet_2_mv2c10_binned"] = static_cast<double>(hasBTag()(leadingsubjets[1]));
      input_vals["sjet_2_e"] = leadingsubjets[1].E()/MeV;
      input_vals["sjet_2_phi"] = leadingsubjets[1].phi();
      input_vals["sjet_2_eta"] = leadingsubjets[1].eta();
      input_vals["sjet_2_pt"] = leadingsubjets[1].pt()/MeV;

      if (leadingsubjets.size() >= 3){
        input_vals["sjet_3_mv2c10_binned"] = static_cast<double>(hasBTag()(leadingsubjets[2]));
        input_vals["sjet_3_e"] = leadingsubjets[2].E()/MeV;
        input_vals["sjet_3_phi"] = leadingsubjets[2].phi();
        input_vals["sjet_3_eta"] = leadingsubjets[2].eta();
        input_vals["sjet_3_pt"] = leadingsubjets[2].pt()/MeV;
      }
      else {
        input_vals["sjet_3_mv2c10_binned"] = -1.0;
        input_vals["sjet_3_e"] = 0.0;
        input_vals["sjet_3_phi"] = totag.phi();
        input_vals["sjet_3_eta"] = totag.eta();
        input_vals["sjet_3_pt"] = 0.0;
      }
    }
    else {
      input_vals["sjet_2_mv2c10_binned"] = -1.0;
      input_vals["sjet_2_e"] = 0.0;
      input_vals["sjet_2_phi"] = totag.phi();
      input_vals["sjet_2_eta"] = totag.eta();
      input_vals["sjet_2_pt"] = 0.0;

      input_vals["sjet_3_mv2c10_binned"] = -1.0;
      input_vals["sjet_3_e"] = 0.0;
      input_vals["sjet_3_phi"] = totag.phi();
      input_vals["sjet_3_eta"] = totag.eta();
      input_vals["sjet_3_pt"] = 0.0;
    }

    compute(input_vals, scoresOut);
  }

  MCBot_TagType MCBot_tagger::tag(const PseudoJet& totag, const Jets &constits){

    std::map<string, double> outputs;
    computeScores(totag, constits, outputs);

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

    MSG_DEBUG("(DV, DH, Dtop, Dlight) = (" << outputs["dnnOutput_V"] << ", " << 
        outputs["dnnOutput_H"] << ", " << outputs["dnnOutput_top"] << ", " <<
         outputs["dnnOutput_light"]  << ")");
    MSG_DEBUG("(PV, PH, Ptop) = (" << PV << ", " << PH << ", " << Ptop << ")");
    MSG_DEBUG("(isV, isH, istop) = (" << isV << ", " << isH << ", " << istop << ")");

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

  void MCBot_tagger::dumpJetToCSV(const PseudoJet& totag, const Jets &constits, 
                                  std::map<string, double>& scoresOut, std::string& filename,
                                  bool overWrite){

    Jets leadingsubjets = sortByPt(constits);                                    
    
    std::ofstream file;
    file.open(filename, overWrite ? std::ofstream::out : std::ofstream::app);
    file << totag.E()*1000 << ", " << totag.m()*1000 << ", " << totag.eta() << ", " <<
            totag.phi() << ", " << leadingsubjets.size() << ", " <<

            leadingsubjets[0].pT()*1000 << ", " << leadingsubjets[0].E()*1000 << ", " <<
            leadingsubjets[0].phi() << ", " << leadingsubjets[0].eta() << ", " <<
            static_cast<double>(hasBTag()(leadingsubjets[0])) << ", ";
    if (leadingsubjets.size() >= 2){
      file << leadingsubjets[1].pT()*1000 << ", " << leadingsubjets[1].E()*1000 << ", " <<
            leadingsubjets[1].phi() << ", " << leadingsubjets[1].eta() << ", " <<
            static_cast<double>(hasBTag()(leadingsubjets[1])) << ", ";

      if (leadingsubjets.size() >= 3){
      file << leadingsubjets[2].pT()*1000 << ", " << leadingsubjets[2].E()*1000 << ", " <<
            leadingsubjets[2].phi() << ", " << leadingsubjets[2].eta() << ", " <<
            static_cast<double>(hasBTag()(leadingsubjets[2])) << ", ";
      } else {
        file << 0.0 << ", " << 0.0 << ", " << totag.phi() << ", " << totag.eta() << ", " <<
            -1.0 << ", ";
      }
    } else {
      file << 0.0 << ", " << 0.0 << ", " << totag.phi() << ", " << totag.eta() << ", " <<
            -1.0 << ", ";
      file << 0.0 << ", " << 0.0 << ", " << totag.phi() << ", " << totag.eta() << ", " <<
            -1.0 << ", ";
    } 
    file << "\n";
    file.close();
  }

  Log& MCBot_tagger::getLog() const {
    return Rivet::Log::getLog("Rivet.MCBot_tagger");
  }





}    