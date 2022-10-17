// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Tools/RivetLWTNN.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ATLAS_2021_I1869040 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2021_I1869040);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9);

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      // TODO: Smear jets
      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");

      // FinalState of direct photons and bare muons and electrons in the event
      DirectFinalState photons(Cuts::abspid == PID::PHOTON);
      DirectFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);

      // Dress the bare direct leptons with direct photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons
      Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
      DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
      declare(dressed_leps, "leptons");

      // Missing momentum
      declare(MissingMomentum(fs), "MET");

      // Initialise the LWTNN graphs (one for each numJets 4-8)
      // Doing this in init saves loading it separately for each event.
      for (size_t i = 4; i < 9; ++i){
        _lwgs[i] = mkGraphLWTNN(analysisDataPath(std::to_string(i)+"j.json"));
      }


      //Book counters and histos.
      //2lsc channel - jet counting
      book(_c["ss_10j20_0b"], "ss_10j20_0b");
      book(_c["ss_10j20_3b"], "ss_10j20_3b");
      book(_c["ss_8j40_0b"], "ss_8j40_0b");
      book(_c["ss_8j40_3b"], "ss_8j40_3b");
      book(_c["ss_7j60_0b"], "ss_7j60_0b");
      book(_c["ss_7j60_3b"], "ss_7j60_3b");
      book(_c["ss_6j80_0b"], "ss_6j80_0b");
      book(_c["ss_6j80_3b"], "ss_6j80_3b");
      book(_c["ss_6j100_0b"], "ss_6j100_0b");
      book(_c["ss_6j100_3b"], "ss_6j100_3b");
      //2lsc channel - shape
      book(_c["ss_shape_6j_3b"], "ss_shape_6j_3b");

      //1l channel - jet counting
      book(_c["1l_15j20_0b"], "1l_15j20_0b");
      book(_c["1l_15j20_3b"], "1l_15j20_3b");
      book(_c["1l_12j40_0b"], "1l_12j40_0b");
      book(_c["1l_12j40_3b"], "1l_12j40_3b");
      book(_c["1l_12j60_0b"], "1l_12j60_0b");
      book(_c["1l_12j60_3b"], "1l_12j60_3b");
      book(_c["1l_10j80_0b"], "1l_10j80_0b");
      book(_c["1l_10j80_3b"], "1l_10j80_3b");
      book(_c["1l_8j100_0b"], "1l_8j100_0b");
      book(_c["1l_8j100_3b"], "1l_8j100_3b");
      
      //1l channel - shape
      //NeuralNetwork Histos:
      book(_h["1l_nn_4j_4b"], "1l_nn_4j_4b", 4, 0, 1);
      book(_h["1l_nn_5j_4b"], "1l_nn_5j_4b", 4, 0, 1);
      book(_h["1l_nn_6j_4b"], "1l_nn_6j_4b", 4, 0, 1);
      book(_h["1l_nn_7j_4b"], "1l_nn_7j_4b", 4, 0, 1);
      book(_h["1l_nn_8j_4b"], "1l_nn_8j_4b", 4, 0, 1);

      //NeuralNetwork sig regions:
      book(_c["1l_shape_4j_4b"], "1l_shape_4j_4b");
      book(_c["1l_shape_5j_4b"], "1l_shape_5j_4b");
      book(_c["1l_shape_6j_4b"], "1l_shape_6j_4b");
      book(_c["1l_shape_7j_4b"], "1l_shape_7j_4b");
      book(_c["1l_shape_8j_4b"], "1l_shape_8j_4b");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /////////////////////////////////////////////////////////////////////////////
      //OBJECT RETRIEVAL
      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 4.5);
      // Retrieve dressed leptons, sorted by pT
      Particles unfiltered_leptons = sortByPt(apply<FinalState>(event, "leptons").particles());
      //Extract electrons from leptons:
      Particles electrons = filter_select(unfiltered_leptons, Cuts::abspid == PID::ELECTRON && Cuts::abseta < 2.47 && Cuts::pt > 15*GeV);
      //Extract muons from leptons:
      Particles muons = filter_select(unfiltered_leptons, Cuts::abspid == PID::MUON && Cuts::abseta < 2.5 && Cuts::pt > 15*GeV);
      ThreeMomentum met = apply<MissingMomentum>(event, "MET").vectorEtMiss();//TODO -> double check MET definitions

      /////////////////////////////////////////////////////////////////////////////
      //OVERLAP REMOVAL (based purely on SimpleAnalaysis script)
      
      //electrons within deltaR 0.01 of a muon
      idiscardIfAnyDeltaRLess(electrons, muons, 0.01);

      //non- btagged jets within 0.2 of an eletron
      idiscardIfAny(jets, electrons,
                      [](const Jet& j, const Particle& e){return ((deltaR(j, e) < 0.4) && !j.bTagged());});
      //TODO -> different b tag criteria here?

      //Remove non b-tagged jets of <3 tracks.
      //The additional check in the SimpleAnalysis for muon.pt > 0.5*jet.pt() is NOT explicitly in the paper.
      //Possible approximation for some extra BDT stuff? Who knows.
      idiscardIfAny(jets, muons, [](const Jet& j,  const Particle &m)
                          {return (!j.bTagged() && (j.size() < 3) && deltaR(j, m) < 0.4);});

      // Electrons and Muons close to jets - Following what paper says
      //electrons within 0.4 of a jet
      idiscardIfAny(electrons, jets, [](const Particle& e, const Jet& j){return deltaR(e, j) < min(0.4, 0.04 + 10*GeV/e.pt());});
      //muons within 0.4 of a jet:
      idiscardIfAny(muons, jets, [](const Particle& m, const Jet& j){return deltaR(m, j) < min(0.4, 0.04 + 10*GeV/m.pt());});
      // Electrons and Muons close to jets - SimpleAnalysis version
      // //electrons within 0.4 of a jet
      // idiscardIfAnyDeltaRLess(electrons, jets, 0.4);
      // //muons within 0.4 of a jet:
      // idiscardIfAnyDeltaRLess(muons, jets, 0.4);
      

      //Derived particle sets after overlap removal
      Particles leptons = sortByPt(electrons + muons);
      // Select jets ghost-associated to B-hadrons with a certain fiducial selection
      Jets bjets = filter_select(jets, hasBTag(Cuts::pT > 5*GeV && Cuts::abseta < 2.5));
      //TODO - x2 check this^


      /////////////////////////////////////////////////////////////////////////////
      //PRESELECTION
      
      //require at least 1 lep, pT >= 27 GeV
      if (leptons.size() == 0 || leptons[0].pt() < 27){
        vetoEvent;
      }
      //require at least four jets, pT > 20
      if (jets.size() < 4){
        vetoEvent;
      }


      ////////////////////////////////////////////////////////////////////////////
      //DECIDE LEPTON CATEGORY

      //TODO - this is as thorough as SimpleAnalysis but not as thorough as the paper.
      if (leptons.size() == 2 && leptons[0].charge()*leptons[1].charge() > 0){
        _is2lsc = true;
      }
      else {
        _is2lsc = false;
      }

      //std::cout << __FILE__ << ": " << __LINE__ << std::endl;

      ////////////////////////////////////////////////////////////////////////////
      //2lsc channel regions:
      if (_is2lsc){
        ////////////////////////////////////////////////////////////////////////////
        //2lsc Jet counting
        const size_t njets40 = std::count_if(jets.begin(), jets.end(),
                                    [](const Jet& j){return (j.pt() >= 40*GeV);});
        const size_t njets60 = std::count_if(jets.begin(), jets.end(), 
                                    [](const Jet& j){return (j.pt() >= 60*GeV);});
        const size_t njets80 = std::count_if(jets.begin(), jets.end(), 
                                    [](const Jet& j){return (j.pt() >= 80*GeV);});
        const size_t njets100 = std::count_if(jets.begin(), jets.end(), 
                                    [](const Jet& j){return (j.pt() >= 100*GeV);});
        const size_t nbjets40 = std::count_if(bjets.begin(), bjets.end(),
                                    [](const Jet& j){return (j.pt() >= 40*GeV);});
        const size_t nbjets60 = std::count_if(bjets.begin(), bjets.end(), 
                                    [](const Jet& j){return (j.pt() >= 60*GeV);});
        const size_t nbjets80 = std::count_if(bjets.begin(), bjets.end(), 
                                    [](const Jet& j){return (j.pt() >= 80*GeV);});
        const size_t nbjets100 = std::count_if(bjets.begin(), bjets.end(), 
                                    [](const Jet& j){return (j.pt() >= 100*GeV);});                                    

        if (jets.size() >= 10){
          if (bjets.size() == 0){
            _c["ss_10j20_0b"]->fill();
          } else if (bjets.size() >= 3){
            _c["ss_10j20_3b"]->fill();
          } 
        }
        if (njets40 >= 8){
          if (nbjets40 == 0){
            _c["ss_8j40_0b"]->fill();
          } else if (nbjets40 >= 3){
            _c["ss_8j40_3b"]->fill();
          } 
        }
        if (njets60 >= 7){
          if (nbjets60 == 0){
            _c["ss_7j60_0b"]->fill();
          } else if (nbjets60 >= 3){
            _c["ss_7j60_3b"]->fill();
          } 
        }
        if (njets80 >= 6){
          if (nbjets80 == 0){
            _c["ss_6j80_0b"]->fill();
          } else if (nbjets80 >= 3){
            _c["ss_6j80_3b"]->fill();
          } 
        }
        if (njets100 >= 6){
          if (nbjets100 == 0){
            _c["ss_6j100_0b"]->fill();
          } else if (nbjets100 >= 3){
            _c["ss_6j100_3b"]->fill();
          } 
        }
        ////////////////////////////////////////////////////////////////////////////
        //2lsc shape analysis
        if (jets.size() == 6 && bjets.size() > 2){
          if (calc_mlj_pair(leptons[0], leptons[1], jets, 4) < 155){
            _c["ss_shape_6j_3b"]->fill();
          }
        }
      }
      ////////////////////////////////////////////////////////////////////////////
      //1l channel
      else {
        ////////////////////////////////////////////////////////////////////////////
        //1l JET COUNTING ANALYSIS
        const size_t njets40 = std::count_if(jets.begin(), jets.end(),
                                    [](const Jet& j){return (j.pt() >= 40*GeV);});
        const size_t njets60 = std::count_if(jets.begin(), jets.end(), 
                                    [](const Jet& j){return (j.pt() >= 60*GeV);});
        const size_t njets80 = std::count_if(jets.begin(), jets.end(), 
                                    [](const Jet& j){return (j.pt() >= 80*GeV);});
        const size_t njets100 = std::count_if(jets.begin(), jets.end(), 
                                    [](const Jet& j){return (j.pt() >= 100*GeV);});
        const size_t nbjets40 = std::count_if(bjets.begin(), bjets.end(),
                                    [](const Jet& j){return (j.pt() >= 40*GeV);});
        const size_t nbjets60 = std::count_if(bjets.begin(), bjets.end(), 
                                    [](const Jet& j){return (j.pt() >= 60*GeV);});
        const size_t nbjets80 = std::count_if(bjets.begin(), bjets.end(), 
                                    [](const Jet& j){return (j.pt() >= 80*GeV);});
        const size_t nbjets100 = std::count_if(bjets.begin(), bjets.end(), 
                                    [](const Jet& j){return (j.pt() >= 100*GeV);});
        if (jets.size() >= 15){
          if (bjets.size() == 0){
            _c["1l_15j20_0b"]->fill();
          } else if (bjets.size() >= 3){
            _c["1l_15j20_3b"]->fill();
          } 
        }
        if (njets40 >= 12){
          if (nbjets40 == 0){
            _c["1l_12j40_0b"]->fill();
          } else if (nbjets40 >= 3){
            _c["1l_12j40_3b"]->fill();
          } 
        }
        if (njets60 >= 12){
          if (nbjets60 == 0){
            _c["1l_12j60_0b"]->fill();
          } else if (nbjets60 >= 3){
            _c["1l_12j60_0b"]->fill();
          } 
        }
        if (njets80 >= 10){
          if (nbjets80 == 0){
            _c["1l_10j80_0b"]->fill();
          } else if (nbjets80 >= 3){
            _c["1l_10j80_0b"]->fill();
          } 
        }
        if (njets100 >= 8){
          if (nbjets100 == 0){
            _c["1l_8j100_0b"]->fill();
          } else if (nbjets100 >= 3){
            _c["1l_8j100_3b"]->fill();
          } 
        }

        ////////////////////////////////////////////////////////////////////////////
        //Neural Network (shape?) based analysis
        if (jets.size() <= 8 && bjets.size() >= 1){
          //I hope this is the only scaling they did as its the only scaling I could find.
          constexpr double Escale = 100 * GeV;

          const map<string, double> nn_input_map = {
            //Whole event level features
            {"n_jet", jets.size()},
            {"n_bcat", (bjets.size() <= 3) ? bjets.size() - 1 : 3},
            {"n_bjet", bjets.size()},
            {"three_jet_maxpt_mass", calc_threejet_max_pt_mass(jets) / Escale},
            {"ht_jet", std::accumulate(jets.begin(), jets.end(), 0.0, 
                          [](double t, const Jet& j){return t + j.pt();}) / Escale},
            {"ht_bjet", std::accumulate(bjets.begin(), bjets.end(), 0.0,
                          [](double t, const Jet& bj){return t + bj.pt();}) / Escale},
            {"three_jet_lep_met_mass_max_pt", 
                calc_threejet_lepmet_max_pt_mass(jets, leptons[0], met)/Escale},
            //TODO see function defintion^^^
            {"minDeltaRlj", min_dr_lep_jet(jets, leptons)},
            {"minmax_mass", calc_minmax_mass(jets) / Escale},
            {"met", met.pt() / Escale},
            {"met_phi", met.phi()},

            //Jet pts (n.b at least four jets must exist to pass presel)
            {"jet_pt_0", jets[0].pt()/Escale},
            {"jet_pt_1", jets[1].pt()/Escale},
            {"jet_pt_2", jets[2].pt()/Escale},
            {"jet_pt_3", jets[3].pt()/Escale},
            {"jet_pt_4", jets.size() > 4 ? jets[4].pt()/Escale: 1e-7},
            {"jet_pt_5", jets.size() > 5 ? jets[5].pt()/Escale: 1e-7},
            {"jet_pt_6", jets.size() > 6 ? jets[6].pt()/Escale: 1e-7},
            {"jet_pt_7", jets.size() > 7 ? jets[7].pt()/Escale: 1e-7},
            {"jet_pt_8", jets.size() > 8 ? jets[8].pt()/Escale: 1e-7},
            {"jet_pt_9", jets.size() > 9 ? jets[9].pt()/Escale: 1e-7},

            //Jet etas (n.b at least four jets must exist to pass presel)
            {"jet_eta_0", jets[0].eta()},
            {"jet_eta_1", jets[1].eta()},
            {"jet_eta_2", jets[2].eta()},
            {"jet_eta_3", jets[3].eta()},
            {"jet_eta_4", jets.size() > 4 ? jets[4].eta(): 1e-7},
            {"jet_eta_5", jets.size() > 5 ? jets[5].eta(): 1e-7},
            {"jet_eta_6", jets.size() > 6 ? jets[6].eta(): 1e-7},
            {"jet_eta_7", jets.size() > 7 ? jets[7].eta(): 1e-7},
            {"jet_eta_8", jets.size() > 8 ? jets[8].eta(): 1e-7},
            {"jet_eta_9", jets.size() > 9 ? jets[9].eta(): 1e-7},

            //Jet phis (n.b at least four jets must exist to pass presel)
            {"jet_phi_0", jets[0].phi()},
            {"jet_phi_1", jets[1].phi()},
            {"jet_phi_2", jets[2].phi()},
            {"jet_phi_3", jets[3].phi()},
            {"jet_phi_4", jets.size() > 4 ? jets[4].phi(): -5},
            {"jet_phi_5", jets.size() > 5 ? jets[5].phi(): -5},
            {"jet_phi_6", jets.size() > 6 ? jets[6].phi(): -5},
            {"jet_phi_7", jets.size() > 7 ? jets[7].phi(): -5},
            {"jet_phi_8", jets.size() > 8 ? jets[8].phi(): -5},
            {"jet_phi_9", jets.size() > 9 ? jets[9].phi(): -5},

            //Jet energies (n.b at least four jets must exist to pass presel)
            {"jet_e_0", jets[0].E()/Escale},
            {"jet_e_1", jets[1].E()/Escale},
            {"jet_e_2", jets[2].E()/Escale},
            {"jet_e_3", jets[3].E()/Escale},
            {"jet_e_4", jets.size() > 4 ? jets[4].E()/Escale: 1e-7},
            {"jet_e_5", jets.size() > 5 ? jets[5].E()/Escale: 1e-7},
            {"jet_e_6", jets.size() > 6 ? jets[6].E()/Escale: 1e-7},
            {"jet_e_7", jets.size() > 7 ? jets[7].E()/Escale: 1e-7},
            {"jet_e_8", jets.size() > 8 ? jets[8].E()/Escale: 1e-7},
            {"jet_e_9", jets.size() > 9 ? jets[9].E()/Escale: 1e-7},

            //Jet btags (n.b at least four jets must exist to pass presel)
            //TODO: WORK OUT
            {"jet_btag_bin_0", (jets[0].bTagged() ? 5 : 1)},
            {"jet_btag_bin_1", (jets[1].bTagged() ? 5 : 1)},
            {"jet_btag_bin_2", (jets[2].bTagged() ? 5 : 1)},
            {"jet_btag_bin_3", (jets[3].bTagged() ? 5 : 1)},
            {"jet_btag_bin_4", jets.size() > 4 ? (jets[4].bTagged() ? 5 : 1): 0},
            {"jet_btag_bin_5", jets.size() > 5 ? (jets[5].bTagged() ? 5 : 1) : 0},
            {"jet_btag_bin_6", jets.size() > 6 ? (jets[6].bTagged() ? 5 : 1) : 0},
            {"jet_btag_bin_7", jets.size() > 7 ? (jets[7].bTagged() ? 5 : 1) : 0},
            {"jet_btag_bin_8", jets.size() > 8 ? (jets[8].bTagged() ? 5 : 1) : 0},
            {"jet_btag_bin_9", jets.size() > 9 ? (jets[9].bTagged() ? 5 : 1) : 0},

            //Leptons features
            {"lepton_pt", leptons[0].pt()},
            {"lepton_eta", leptons[0].eta()},
            {"lepton_phi", leptons[0].phi()},
            {"lepton_e", leptons[0].E()},
          };
          map<string, map<string, double>> nn_input = {{"node_0", nn_input_map}};
          //Compute score using netwrok corresponding to current numJets.
          map<string, double> nn_output = _lwgs[jets.size()]->compute(nn_input);

          //Fill appropriate histo:
          _h["1l_nn_"+std::to_string(jets.size())+"j_4b"] -> fill(nn_output["out_0"]);
  
          //If >= 4 bjets, and NN cut achieved, fill appropriate sigregion also:
          //If >= 4 bjets, and NN cut achieved, fill appropriate sigregion also:
          if(nn_output["out_0"] > _nnCuts[jets.size()])
            _c["1l_shape_"+std::to_string(jets.size())+"j_4b"]->fill();

        }
      }



      //std::cout << __FILE__ << ": " << __LINE__ << std::endl;

    }


    /// Normalise histograms etc., after the run
    void finalize() {

    }

    /// @}


    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, CounterPtr> _c;
    /// @}

    

    private:
    /// LWTNN neural-network pointers (one for each case 4jets - 8jets)
    map<size_t, unique_ptr<lwt::LightweightGraph>> _lwgs;

    ///Store whether event belongs in 1l or 2lsc region:
    bool _is2lsc;

    ///Cut-off values for the neural network, as a map of {njets, cut}
    map<size_t, double> _nnCuts = {{4, 0.73}, {5, 0.76}, {6, 0.77}, {7, 0.72}, {8, 0.73}};
    

    ///@name utility functions for this analysis
    /// @{
    //Get the mass of the system of the three jets with highest combined pt
    double calc_threejet_max_pt_mass(const Jets& jets){
      double max_pt = 0;
      FourMomentum max4mom;
      for (size_t i = 0; i < jets.size(); ++i){
        for (size_t j = i+1; j < jets.size(); ++j){
          for (size_t k = j+1; k < jets.size(); ++k){
            if ((jets[i].mom() + jets[j].mom() + jets[k].mom()).pt() > max_pt){
              max_pt = (jets[i].mom() + jets[j].mom() + jets[k].mom()).pt();
              max4mom = (jets[i].mom() + jets[j].mom() + jets[k].mom());
            }
          } 
        }  
      }
      return max4mom.mass();
    }

    //Get the mass of the system of the three jets with highest combined pt
    //TODO: WORK OUT WTF is going on with MET -> ignored for now.
    double calc_threejet_lepmet_max_pt_mass(const Jets& jets, const Particle &lep,
                                            const ThreeMomentum met){
      double max_pt = 0;
      FourMomentum max4mom;
      for (size_t i = 0; i < jets.size(); ++i){
        for (size_t j = i+1; j < jets.size(); ++j){
          for (size_t k = j+1; k < jets.size(); ++k){
            if (((ThreeMomentum)jets[i].p3() + (ThreeMomentum)jets[j].p3() 
                    + (ThreeMomentum)jets[k].p3() + (ThreeMomentum)lep.p3()
                      + met).pt() > max_pt){
              max_pt = (jets[i].mom() + jets[j].mom() + jets[k].mom() + lep.mom()).pt();
              max4mom = (jets[i].mom() + jets[j].mom() + jets[k].mom() + lep.mom());
            }
          } 
        }  
      }
      return max4mom.mass();
    }

    //Calculate the minium deltaR
    double min_dr_lep_jet(const Jets& jets, const Particles& leptons){
      double min_dr = DBL_MAX;
      for (const Particle &l : leptons){
        for (const Jet &j: jets){
          if (deltaR(j, l) < min_dr){
            min_dr = deltaR(j, l);
          }
        }
      }
      return min_dr;
    }

    //Helper function for the below, copied from SA for same reason.
    int countSetBits( int n)
    {
      int count = 0;
      while (n) {
          count += n & 1;
          n >>= 1;
      }
      return count;
    }

    //TODO (or at least nb) Copied almost wholesale from simpleanalysis as 
    //I cannot understand thi
    double calc_minmax_mass(const Jets& jets, int jetdiff=10)
    {
      const int nJets = jets.size();

      //bitwise representation of all jets and for which half they are selected
      // One bit for every jet, marking into which set they are grouped
      const unsigned int bitmax = 1 << nJets;
      float minmass = 999999999;

      for(unsigned int bit=0; bit < bitmax; bit++){
          const int bitcount = countSetBits(bit);
          if (abs(nJets - 2*bitcount) > jetdiff) {
            continue;
          }

          FourMomentum sum1, sum2;
          // loop through jets and assign to either sum, depending on bit
          for(int i=0; i<nJets; i++) {
              if (bit & (1<<i)) 
                sum1 += jets[i];
              else
                sum2 += jets[i];
          }
          if (sum1.mass() > sum2.mass() && sum1.mass() < minmass) 
            minmass = sum1.mass();
          else if (sum2.mass() > sum1.mass() && sum2.mass() < minmass)
            minmass = sum2.mass();
      }

      return minmass;
    }


    //Copied wholesale from SimpleAnalysis with some logic tweaks
    static double calc_mlj_pair(const Particle& l1, const Particle& l2,
                            const Jets& jets, const size_t max_njets)
    {
      float minmax = 9e9;
      float tmpmass;
      const size_t iterate_max = max(max_njets, jets.size());
      for(size_t i = 0; i < iterate_max; ++i){
        for(size_t j = 0; j < iterate_max; ++j){
          if(i==j) continue;
          tmpmass = std::max((l1.momentum() + jets[i].momentum()).mass(),
                              (l2.momentum() + jets[j].momentum()).mass());
          if(tmpmass < minmax) minmax = tmpmass;
        }
      }
      return minmax;
    }

    

    /// @}

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2021_I1869040);

}
