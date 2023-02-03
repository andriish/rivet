// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Tools/RivetLWTNN.hh"
#include "Rivet/Projections/SmearedParticles.hh"
#include "Rivet/Projections/SmearedJets.hh"

namespace Rivet {

  enum MCBot_TagType{
    Bkg,
    Vec,
    Higgs,
    top
  };

  /// @brief Add a short analysis description here
  class ATLAS_2022_I2172216 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2022_I2172216);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      _mode = 1;
      const std::string mode = getOption("MODE");
      if (mode == "VAL-SIG"){
        _mode = 2;
      }
      else if (mode == "VAL-BKG"){
        _mode = 3;
      }


      // Initialise and register projections
    
      //Final state of all tracks - used in overlap removal
      declare(VisibleFinalState(Cuts::abseta < 4.9), "allTracks");

      // electrons
      // n.b. cut too low to accounnt for smearing
      const DirectFinalState bare_leptons((Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON));
      const DirectFinalState photons(Cuts::abspid == PID::PHOTON);
      DressedLeptons dressed_leps(photons, bare_leptons, 0.1, Cuts::abseta < 2.47 && Cuts::pt > 20*GeV && !Cuts::absetaIn(1.37, 1.52));

	    SmearedParticles recoelectrons(dressed_leps, ELECTRON_EFF_ATLAS_RUN2_TIGHT, ELECTRON_SMEAR_ATLAS_RUN2, Cuts::abspid==PID::ELECTRON);
    	declare(recoelectrons, "SmearedElec");
      /// muons
      // n.b. cut too low to account for smearing - will cut again in preselection
	    SmearedParticles recomuons(dressed_leps, MUON_EFF_ATLAS_RUN2, MUON_SMEAR_ATLAS_RUN2, Cuts::abspid==PID::MUON);
    	declare(recomuons, "SmearedMuon");

      /// jets
	    const FinalState fsj(Cuts::abseta < 4.8);
	    FastJets Sj(fsj, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE);
	    declare("Sjet", Sj);
	    SmearedJets SSj(Sj, JET_SMEAR_ATLAS_RUN2, JET_BTAG_PERFECT);
      /// @todo Also look into angular smearing? Need a custom smearing function, building on the ATLAS R2
      //SmearedJets SSj(Sj, JET_SMEAR_ANGULAR, JET_BTAG_EFFS(0.77, 1./6.2, 1./134));
    	declare(SSj, "smearedSjet");

      declare(MissingMom(), "Etmiss");

      //Load the NN:
      _nn = mkLWTNN("/home/tomek/PHYSICS_INSTALLS/rivet_seven/rivet/analyses/examples/ATLAS_2018_I1685207.nn.json");


      // Book histograms and counters.
      //TODO: This will be a lot easier when they finally make a hepdata entry
      //(IF all auxiliary plots are included!)

      /// VALIDATION:
      if (_mode > 1){
        book(_h["DNN_V"], "DNN_V", 30, 0., 1.);
        book(_h["DNN_H"], "DNN_H", 30, 0., 1.);
        book(_h["DNN_top"], "DNN_top", 30, 0., 1.);
        
        if (_mode == 2){
          book(_h["DNN_V_forH"], "DNN_V_forH", 30, 0., 1.);
          book(_h["DNN_V_fortop"], "DNN_V_fortop", 30, 0., 1.);
          
          book(_h["DNN_H_forV"], "DNN_H_forV", 30, 0., 1.);
          book(_h["DNN_H_fortop"], "DNN_H_fortop", 30, 0., 1.);

          book(_h["DNN_top_forV"], "DNN_top_forV", 30, 0., 1.);
          book(_h["DNN_top_forH"], "DNN_top_forH", 30, 0., 1.);

          book(_h["DNN_light_forV"], "DNN_light_forV", 30, 0., 1.);
          book(_h["DNN_light_forH"], "DNN_light_forH", 30, 0., 1.);
          book(_h["DNN_light_fortop"], "DNN_light_fortop", 30, 0., 1.);
        }
        else if (_mode == 3){
          book(_h["DNN_light"], "DNN_light", 30, 0., 1.);
        }
        

        
        return;
      }

      /// SIGNAL:
      //2L 1b SR:
      for (const string& s : _2L_MCBOT_categories){
        book(_h["Hist2l_1b_SR_"+s], "Hist2l_1b_SR_"+s, {0, 600, 1000, 1400, 2000} );
        book(_sigBins["2l_1b_SR_"+s], "2l_1b_SR_"+s);
      }
      // 2l 2b SR:
      for (const string& s : _2L_MCBOT_categories){
        book(_h["Hist2l_2b_SR_"+s], "Hist2l_2b_SR_"+s, 
          (s == "Notag" || s =="Htag" || s=="toptag") ? 
            // I was today years old when I learned you can't put an initialiser list
            // directly inside a ternary operator.
            (vector<double>{0, 600, 1000, 1400, 2000}) : vector<double>{0, 600, 1000, 2000} );
        book(_sigBins["2l_2b_SR_"+s], "2l_2b_SR_"+s);
      }
      // 3l SR
      // The fact that there are so many binnings here shows how useful hepdata would be!
      book(_h["Hist3l_SR_Htag"], "Hist3l_SR_Htag", {1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400});
      book(_h["Hist3l_SR_OF"], "Hist3l_SR_OF", {1000, 1200, 1400, 1600, 1800, 2000, 2400});
      book(_h["Hist3l_SR_Notag"], "Hist3l_SR_Notag", {1000, 1200,  2400});
      book(_h["Hist3l_SR_toptag"], "Hist3l_SR_toptag", {1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400});
      book(_h["Hist3l_SR_Vtag"], "Hist3l_SR_Vtag", {1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400});
      for (const string& s : _3L_MCBOT_categories){
        book(_sigBins["3l_SR_"+s], "3l_SR_"+s);
      }

      /// CONTROL:
      book(_h["Hist2l_1b_CR"], "Hist2l_1b_CR", {920, 1150, 1380});
      book(_controlBins["2l_1b_CR"], "2l_1b_CR");
      book(_h["Hist2l_2b_CR"], "Hist2l_2b_CR", {920, 1150, 1380});
      book(_controlBins["2l_2b_CR"], "2l_2b_CR");
      book(_h["Hist3l_VV_CR"], "Hist3l_VV_CR", 16, 300, 1500);
      book(_controlBins["3l_VV_CR"], "3l_VV_CR");

      /// DEBUG:
      book(_h["2lEtmiss"], "2lEtmiss", 40, 0, 4000);
      book(_h["2lHTjet"], "2lHTjet",  40, 0, 4000);
      book(_h["2lHTjetcentral"],"2lHTjetcentral", 40, 0, 4000);
      book(_h["2lEtmiss_HTjet"],"2lEtmiss_HTjet", 40, 0, 4000);
      book(_h["2lEtmiss_HTjetcentral"], "2lEtmiss_HTjetcentral", 40, 0, 4000);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      /// OBJECTRETRIEVAL:
      // Retrieve Leptons
      Particles smeared_electrons = apply<ParticleFinder>(event, "SmearedElec").particles(
        Cuts::abseta < 2.47 && !Cuts::absetaIn(1.37, 1.52) && Cuts::pT > 28*GeV);
      Particles smeared_muons = apply<ParticleFinder>(event, "SmearedMuon").particles(
        Cuts::abseta < 2.5 && Cuts::pt > 28*GeV);
      Particles alltracks = apply<ParticleFinder>(event, "allTracks").particles(); //Needed for muon overlap removal
      //Muon overlap removal (e FixedCutTightTrackOnly working point):
      //TODO: The lambdas have gone too far, haven't they...
      idiscard(smeared_muons, [&alltracks]
        (const Particle& mu){
          const double dR = min({0.3, 10*GeV/mu.pt()});
          return (accumulate(alltracks, 0., [&mu, dR](double tot, const Particle &p){
            return deltaR(p, mu) < dR ? tot + p.pt() : tot;
          }) > 1.15*mu.pt());
        });
      
      //Retrieve Jets
      Jets smeared_small_jets = apply<JetAlg>(event, "smearedSjet").jetsByPt(
        (Cuts::pT > 25*GeV && Cuts::abseta < 2.5) || (Cuts::pT > 35*GeV && Cuts::abseta < 4.5)
      );
      //B-tagging -> Note this is a little ugly as we want all 4 working points consistently.
      vector<pair<ThreeMomentum,int>> jet_mv2c10_bands;
      jet_mv2c10_bands.reserve(smeared_small_jets.size());
      for (Jet &j : smeared_small_jets){
        jet_mv2c10_bands.push_back({j.p3(), dummy_mv2c10_band(j)});
      }
      const Jets smeared_bjets = filter_select(smeared_small_jets, hasBTag());
      // ETMiss
      const MissingMomentum ETMiss = apply<MissingMomentum>(event, "Etmiss");

      /// OVERLAPREMOVAL:
      // 1. + 2.  Any muon leaving energy deposits in the calo
      // and "sharing a track" with an electron -> Let's call this dR < 0.05 (TODO: x2 check)
      // Then any electrons sharing track with remaining muons.
      // I don't think this will really happen in MC.
      idiscardIfAnyDeltaRLess(smeared_muons, smeared_electrons, 0.05);
      // 3. Any jet within deltaR 0.2 of an e-
      idiscardIfAnyDeltaRLess(smeared_small_jets, smeared_electrons, 0.2);
      // 4. Any electron within deltaR 0.4 of remaining jets
      idiscardIfAnyDeltaRLess(smeared_electrons, smeared_small_jets, 0.2);
      // 5. Any not b-tagged jet, =< 2 tracks, deltaR 0.2 pf a muon
      idiscardIfAny(smeared_small_jets, smeared_muons,
        [](const Jet &j, const Particle &mu){
          return !j.bTagged() && j.size() < 3 && deltaR(mu.p3(), j.p3()) < 0.2;
        });
      // 6. Any muons within deltaR < min(0.4, 0.04+10/mu.pt) of remaining jet
      idiscardIfAny(smeared_muons, smeared_small_jets,
        [](const Particle &mu, const Jet &j){
          return (deltaR(j.p3(), mu.p3()) < min({0.4, 0.04+10*GeV/mu.pt()}));
        });
      
      /// GET AND TRIM LARGE-R JETS.
      fastjet::JetDefinition jdef(fastjet::antikt_algorithm, 1.0);
      ClusterSequence cseq(smeared_small_jets, jdef);
      PseudoJets large_jets = cseq.inclusive_jets();

      // Trimming
      //Note this is super ugly because we need to keep track of jet substructure by hand.
      PseudoJets TrimmedLargejets;
      vector<PseudoJets> NewConstits;
      naive_trimmer(large_jets, TrimmedLargejets, 0.05, NewConstits);
      PseudoJet tr_pj;

      PseudoJets UnsortedFilteredLargejets;
      vector<PseudoJets> UnsortedFilteredNewConstits;
      for (size_t i = 0; i < TrimmedLargejets.size(); ++i){
        tr_pj = TrimmedLargejets[i];
        if (tr_pj.pt() > 150*GeV && abs(tr_pj.eta()) < 2.5 && tr_pj.m() > 40*GeV){
          UnsortedFilteredLargejets.push_back(tr_pj);
          UnsortedFilteredNewConstits.push_back(NewConstits[i]);
        }
      }

      //Ensure Jets are still pT ordered after trimming:
      // Note that the Consituents need to go with it, which leads to some ugliness
      // Feels like there should be an STL algroithm to make this less ugly.
      PseudoJets FilteredLargejets(UnsortedFilteredLargejets.size());
      vector<PseudoJets> FilteredNewConstits(UnsortedFilteredLargejets.size());
      {
        //Build a permutation vector to get the order we want.
        vector<size_t> permutation(UnsortedFilteredLargejets.size());
        std::iota(permutation.begin(), permutation.end(), 0);
        sort(permutation.begin(), permutation.end(),
         [&UnsortedFilteredLargejets](size_t i, size_t j){
            return UnsortedFilteredLargejets[i].pt() > UnsortedFilteredLargejets[j].pt();
        });
        //Note there's a transform in Rivet namespace also.
        std::transform(permutation.begin(), permutation.end(), 
          FilteredLargejets.begin(), 
            [&UnsortedFilteredLargejets](size_t i){return UnsortedFilteredLargejets[i];
          });
        std::transform(permutation.begin(), permutation.end(), 
          FilteredNewConstits.begin(), [&UnsortedFilteredNewConstits](size_t i){return UnsortedFilteredNewConstits[i];});
      }

      /// VALIDATION:
      if (_mode > 1){
        Particles Vectors;
        Particles Higgses;
        Particles Tops;
        getVHandtop_fromEvent(Vectors, event, {23,24});
        getVHandtop_fromEvent(Higgses, event, {25});
        getVHandtop_fromEvent(Tops, event, {6});

        PseudoJets signal = FilteredLargejets;
        vector<PseudoJets> SignalConstits = FilteredNewConstits;
        vector<MCBot_TagType> JetTags;
        for (size_t counter = 0; counter < signal.size(); ++counter){
          auto& j = signal[counter];
          //First we need a vector of all the constituents, that has b-tagging info. 
          Jets tagged_constituents;
          for (const PseudoJet& pj : SignalConstits[counter]){
            auto it = std::find_if(smeared_small_jets.begin(), smeared_small_jets.end(),
                                    [&pj](Jet& j){return (deltaR(momentum3(pj), j) < 0.1);});
            if (it != smeared_small_jets.end()){
              tagged_constituents.push_back(*it);
            }
            else {
              MSG_WARNING("FAILED TO FIND CONSTITUENT");
            }
          }
          //We need to make sure that the tagged constituents are still pT ordered
          isortByPt(tagged_constituents);

          const size_t nConstits = tagged_constituents.size();

          //Define the MCBot TagType (to be done either by the NN or by efficiencies.
          std::map<string, double> outputs;
          const std::map<string, double> NN_Input = {
            {"rcjet_pt", j.pt()/MeV},
            {"rcjet_numConstituents", static_cast<double>(tagged_constituents.size())},
            {"rcjet_m", j.m()/MeV},
            
            //Lead Jet
            {"sjet_1_mv2c10_binned", get_mv2c10(tagged_constituents[0], jet_mv2c10_bands)},
            {"sjet_1_e", tagged_constituents[0].E()/MeV},
            {"sjet_1_phi", tagged_constituents[0].phi()},
            {"sjet_1_eta", tagged_constituents[0].eta()},
            {"sjet_1_pt", tagged_constituents[0].pt()/MeV},

            //Second Jet
            {"sjet_2_mv2c10_binned", (nConstits > 1) ? get_mv2c10(tagged_constituents[1], jet_mv2c10_bands) : -1. },
            {"sjet_2_e", (nConstits > 1) ? tagged_constituents[1].E()/MeV : 0. },
            {"sjet_2_phi", (nConstits > 1) ? tagged_constituents[1].phi() : j.phi() },
            {"sjet_2_eta", (nConstits > 1) ? tagged_constituents[1].eta() : j.eta() },
            {"sjet_2_pt", (nConstits > 1) ? tagged_constituents[1].pt()/MeV : 0. },

            //Third Jet
            {"sjet_3_mv2c10_binned", (nConstits > 2) ? get_mv2c10(tagged_constituents[2], jet_mv2c10_bands) : -1. },
            {"sjet_3_e", (nConstits > 2) ? tagged_constituents[2].E()/MeV : 0. },
            {"sjet_3_phi", (nConstits > 2) ? tagged_constituents[2].phi() : j.phi() },
            {"sjet_3_eta", (nConstits > 2) ? tagged_constituents[2].eta() : j.eta() },
            {"sjet_3_pt", (nConstits > 2) ? tagged_constituents[2].pt()/MeV : 0. }
          };

          outputs = _nn->compute(NN_Input);
          
          //Signal validation mode
          if (_mode == 2){
            // Work out what signal it is.
            auto itVec = std::find_if(Vectors.begin(), Vectors.end(),
            //TODO: is there a preferred syntax for writing long ugly lambdas?
                        [&j](const Particle& VHTop){
                            return (deltaR(momentum3(VHTop.pseudojet()), momentum3(j)) < 0.75);
                          });
            auto itVec2 = (itVec == Vectors.end())?(Vectors.end()):std::find_if(itVec+1, Vectors.end(),
                        [&j](const Particle& VHTop){
                            return (deltaR(momentum3(VHTop.pseudojet()), momentum3(j)) < 0.75);
                          });

            auto itHiggs = std::find_if(Higgses.begin(), Higgses.end(),
            //TODO: is there a preferred syntax for writing long ugly lambdas?
                        [&j](const Particle& VHTop){
                            return (deltaR(momentum3(VHTop.pseudojet()), momentum3(j)) < 0.75);
                          });
            auto itHiggs2 = (itHiggs == Higgses.end())?(Higgses.end()):std::find_if(itHiggs+1, Higgses.end(),
                        [&j](const Particle& VHTop){
                            return (deltaR(momentum3(VHTop.pseudojet()), momentum3(j)) < 0.75);
                          });

            auto itTop = std::find_if(Tops.begin(), Tops.end(),
            //TODO: is there a preferred syntax for writing long ugly lambdas?
                        [&j](const Particle& VHTop){
                            return (deltaR(momentum3(VHTop.pseudojet()), momentum3(j)) < 0.75);
                          });
            auto itTop2 = (itTop == Tops.end())?(Tops.end()):std::find_if(itTop+1, Tops.end(),
                        [&j](const Particle& VHTop){
                            return (deltaR(momentum3(VHTop.pseudojet()), momentum3(j)) < 0.75);
                          });

            //MSG_INFO((itVec == Vectors.end()) << ", " << (itVec2 == Vectors.end())  << ", " <<  (itHiggs == Higgses.end()) << ", " <<  (itHiggs2 == Higgses.end()) << ", " <<  (itTop == Tops.end()) << ", " <<  (itTop2 == Tops.end()));


            if (itVec != Vectors.end() && itVec2 == Vectors.end() && itHiggs == Higgses.end() && itTop == Tops.end()){
              // This jet matches to a vector:
              _h["DNN_V"]->fill(outputs.at("dnnOutput_V"));
              _h["DNN_H_forV"]->fill(outputs.at("dnnOutput_H"));
              _h["DNN_top_forV"]->fill(outputs.at("dnnOutput_top"));
              _h["DNN_light_forV"]->fill(outputs.at("dnnOutput_light"));
            }
            else if (itVec == Vectors.end() && itHiggs != Higgses.end() && itHiggs2 == Higgses.end() && itTop == Tops.end()){
              // This jet matches to a Higgs:
              _h["DNN_H"]->fill(outputs.at("dnnOutput_H"));
              _h["DNN_V_forH"]->fill(outputs.at("dnnOutput_V"));
              _h["DNN_top_forH"]->fill(outputs.at("dnnOutput_top"));
              _h["DNN_light_forH"]->fill(outputs.at("dnnOutput_light"));
            }
            else if (itVec == Vectors.end() && itHiggs == Higgses.end() && itTop != Tops.end() && itTop2 == Tops.end()){
              // This jet matches to a top:
              _h["DNN_top"]->fill(outputs.at("dnnOutput_top"));
              _h["DNN_V_fortop"]->fill(outputs.at("dnnOutput_V"));
              _h["DNN_H_fortop"]->fill(outputs.at("dnnOutput_H"));
              _h["DNN_light_fortop"]->fill(outputs.at("dnnOutput_light"));
            }
          }
          //Background validation mode:
          else if (_mode == 3){
            _h["DNN_V"]->fill(outputs.at("dnnOutput_V"));
            _h["DNN_H"]->fill(outputs.at("dnnOutput_H"));
            _h["DNN_top"]->fill(outputs.at("dnnOutput_top"));
            _h["DNN_light"]->fill(outputs.at("dnnOutput_light"));
          }
        }
        return;
      }
      
      /// PRESELECTION:

      // 1. 2 opposite sign, same flavour leptons:
      // Get all OSSF pairs:
      vector<pair<Particle, Particle>> OSSF_PAIRS;
      for (size_t i = 0; i < smeared_electrons.size(); ++i){
        for (size_t j = 0; j < i; ++j){
          if (smeared_electrons[i].charge() * smeared_electrons[j].charge() < 1){
            OSSF_PAIRS.push_back({smeared_electrons[i], smeared_electrons[j]});
          }
        }
      }
      for (size_t i = 0; i < smeared_muons.size(); ++i){
        for (size_t j = 0; j < i; ++j){
          if (smeared_muons[i].charge() * smeared_muons[j].charge() < 1){
            OSSF_PAIRS.push_back({smeared_muons[i], smeared_muons[j]});
          }
        }
      }
      if (OSSF_PAIRS.size() < 1)
        vetoEvent;

      // 2. two jets (presumed small?) in central region
      if (count_if(smeared_small_jets.begin(), smeared_small_jets.end(),
        [](const Jet &j){return (j.abseta() < 2.5 );}) < 2){
          vetoEvent;
      }

      // 3. Z-boson candidate |M(ll)-MZ| < 10*GeV
      // Get the Z bsoson candidate - i.e. the OSSF pair with lowest dM
      double best_deltaM = DBL_MAX;
      FourMomentum ZCandidate;
      for (const pair<Particle, Particle>& ossf_pair : OSSF_PAIRS){
        const FourMomentum this_candidate(ossf_pair.first.mom() + ossf_pair.second.mom());
        if (abs(this_candidate.mass()-91.2*GeV) < best_deltaM){
          best_deltaM = abs(this_candidate.mass()-91.2*GeV);
          ZCandidate = this_candidate;
        }
      }
      if (best_deltaM >= 10*GeV)
        vetoEvent;

      /// CHANNELS:

      /// TWOLCHANNEL:
      if (smeared_muons.size() + smeared_electrons.size() == 2){
        /// MORE PRESEL SPECIFIC TO THIS CHANNEL:
        // Z Candidate pT > 300*Gev:
        if (ZCandidate.pt() <= 300)
          vetoEvent;
        //HTjet + ETMiss > 920*GeV
        //TODO: Should this be central only or central + outer?
        const double HTjet = accumulate(smeared_small_jets, 0.0, 
          [](const double tot, const Jet& j){return tot+j.pt();});
        const double HTjetCentral = accumulate(smeared_small_jets, 0.0, 
          [](const double tot, const Jet& j){return j.abseta() < 2.5 ? tot+j.pt() : tot;});          


        _h["2lEtmiss"]->fill(ETMiss.scalarPtMiss());
        _h["2lHTjet"]->fill(HTjet);
        _h["2lHTjetcentral"]->fill(HTjetCentral);
        _h["2lEtmiss_HTjet"]->fill(HTjet + ETMiss.scalarPtMiss());
        _h["2lEtmiss_HTjetcentral"]->fill(HTjetCentral + ETMiss.scalarPtMiss());
        if (HTjet + ETMiss.scalarPtMiss() <= 920*GeV )
          vetoEvent;

        if (smeared_bjets.size() == 0)
            vetoEvent;

        //Populate the control region:
        if (HTjet + ETMiss.scalarPtMiss() <= 1380*GeV ){
          if (smeared_bjets.size() == 1){
            _controlBins["2l_1b_CR"]->fill();
            _h["Hist2l_1b_CR"]->fill(HTjet + ETMiss.scalarPtMiss());
          }
          else {
            _controlBins["2l_2b_CR"]->fill();
            _h["Hist2l_2b_CR"]->fill(HTjet + ETMiss.scalarPtMiss());
          }
          //And now we've done all we need to in the CR.
          return;
        }
        //If not in a control Region, need to run MCBOT tagging.
        PseudoJets signal = FilteredLargejets;
        vector<PseudoJets> SignalConstits = FilteredNewConstits;
        vector<MCBot_TagType> JetTags;
        for (size_t counter = 0; counter < signal.size(); ++counter){
          auto& j = signal[counter];
          //First we need a vector of all the constituents, that has b-tagging info. 
          Jets tagged_constituents;
          for (const PseudoJet& pj : SignalConstits[counter]){
            auto it = std::find_if(smeared_small_jets.begin(), smeared_small_jets.end(),
                                    [&pj](Jet& j){return (deltaR(momentum3(pj), j) < 0.1);});
            if (it != smeared_small_jets.end()){
              tagged_constituents.push_back(*it);
            }
            else {
              MSG_WARNING("FAILED TO FIND CONSTITUENT");
            }
          }
          //We need to make sure that the tagged constituents are still pT ordered
          isortByPt(tagged_constituents);

          const size_t nConstits = tagged_constituents.size();

          //Define the MCBot TagType (to be done either by the NN or by efficiencies.
          MCBot_TagType tag;
          std::map<string, double> outputs;
          const std::map<string, double> NN_Input = {
            {"rcjet_pt", j.pt()/MeV},
            {"rcjet_numConstituents", static_cast<double>(tagged_constituents.size())},
            {"rcjet_m", j.m()/MeV},
            
            //Lead Jet
            {"sjet_1_mv2c10_binned", get_mv2c10(tagged_constituents[0], jet_mv2c10_bands)},
            {"sjet_1_e", tagged_constituents[0].E()/MeV},
            {"sjet_1_phi", tagged_constituents[0].phi()},
            {"sjet_1_eta", tagged_constituents[0].eta()},
            {"sjet_1_pt", tagged_constituents[0].pt()/MeV},

            //Second Jet
            {"sjet_2_mv2c10_binned", (nConstits > 1) ? get_mv2c10(tagged_constituents[1], jet_mv2c10_bands) : -1. },
            {"sjet_2_e", (nConstits > 1) ? tagged_constituents[1].E()/MeV : 0. },
            {"sjet_2_phi", (nConstits > 1) ? tagged_constituents[1].phi() : j.phi() },
            {"sjet_2_eta", (nConstits > 1) ? tagged_constituents[1].eta() : j.eta() },
            {"sjet_2_pt", (nConstits > 1) ? tagged_constituents[1].pt()/MeV : 0. },

            //Third Jet
            {"sjet_3_mv2c10_binned", (nConstits > 2) ? get_mv2c10(tagged_constituents[2], jet_mv2c10_bands) : -1. },
            {"sjet_3_e", (nConstits > 2) ? tagged_constituents[2].E()/MeV : 0. },
            {"sjet_3_phi", (nConstits > 2) ? tagged_constituents[2].phi() : j.phi() },
            {"sjet_3_eta", (nConstits > 2) ? tagged_constituents[2].eta() : j.eta() },
            {"sjet_3_pt", (nConstits > 2) ? tagged_constituents[2].pt()/MeV : 0. }
          };

          outputs = _nn->compute(NN_Input);
          tag =  getTag(outputs);
          JetTags.push_back(tag);
        }
        const string category = get_2l_MCBOT_category(JetTags, smeared_bjets.size());
        _sigBins[(string)"2l_"+((smeared_bjets.size() > 1)?"2":"1")+"b_SR_"+category]->fill(); 
        if (smeared_bjets.size() == 1)       
          _h[(string)"Hist2l_"+((smeared_bjets.size() > 1)?"2":"1")+"b_SR_"+category]->fill((ZCandidate+smeared_bjets[0].mom()).mass());
        else
          _h[(string)"Hist2l_"+((smeared_bjets.size() > 1)?"2":"1")+"b_SR_"+category]->fill((ZCandidate+smeared_bjets[1].mom()).mass());
      } 
      /// THREELCHANNEL:
      else {
        if (ZCandidate.pt() <= 200)
          vetoEvent;
        //HTjetlep > 300*GeV
        //TODO: Should this be central only or central + outer?
        const double HTjetlep = accumulate(smeared_small_jets, 0.0, 
          [](const double tot, const Jet& j){return tot+j.pt();}) + 
          accumulate(smeared_electrons, 0.0, [](const double tot, const Particle& e){
            return tot+e.pt();}) +
          accumulate(smeared_muons, 0.0, [](const double tot, const Particle& mu){
            return tot+mu.pt();});
        if (HTjetlep < 300*GeV )
          vetoEvent;
        
        if (smeared_bjets.size() == 0){
          //We're in the VV_CR
          _controlBins["3l_VV_CR"]->fill();
          _h["Hist3l_VV_CR"]->fill(HTjetlep);
        }
        else {
          // We need to do MCBot tags
          PseudoJets signal = FilteredLargejets;
          vector<PseudoJets> SignalConstits = FilteredNewConstits;
          vector<MCBot_TagType> JetTags;
          for (size_t counter = 0; counter < signal.size(); ++counter){
            auto& j = signal[counter];
            //First we need a vector of all the constituents, that has b-tagging info. 
            Jets tagged_constituents;
            for (const PseudoJet& pj : SignalConstits[counter]){
              auto it = std::find_if(smeared_small_jets.begin(), smeared_small_jets.end(),
                                      [&pj](Jet& j){return (deltaR(momentum3(pj), j) < 0.1);});
              if (it != smeared_small_jets.end()){
                tagged_constituents.push_back(*it);
              }
              else {
                MSG_WARNING("FAILED TO FIND CONSTITUENT");
              }
            }
            //We need to make sure that the tagged constituents are still pT ordered
            isortByPt(tagged_constituents);

            const size_t nConstits = tagged_constituents.size();

            //Define the MCBot TagType (to be done either by the NN or by efficiencies.
            MCBot_TagType tag;
            std::map<string, double> outputs;
            const std::map<string, double> NN_Input = {
              {"rcjet_pt", j.pt()/MeV},
              {"rcjet_numConstituents", static_cast<double>(tagged_constituents.size())},
              {"rcjet_m", j.m()/MeV},
              
              //Lead Jet
              {"sjet_1_mv2c10_binned", get_mv2c10(tagged_constituents[0], jet_mv2c10_bands)},
              {"sjet_1_e", tagged_constituents[0].E()/MeV},
              {"sjet_1_phi", tagged_constituents[0].phi()},
              {"sjet_1_eta", tagged_constituents[0].eta()},
              {"sjet_1_pt", tagged_constituents[0].pt()/MeV},

              //Second Jet
              {"sjet_2_mv2c10_binned", (nConstits > 1) ? get_mv2c10(tagged_constituents[1], jet_mv2c10_bands) : -1. },
              {"sjet_2_e", (nConstits > 1) ? tagged_constituents[1].E()/MeV : 0. },
              {"sjet_2_phi", (nConstits > 1) ? tagged_constituents[1].phi() : j.phi() },
              {"sjet_2_eta", (nConstits > 1) ? tagged_constituents[1].eta() : j.eta() },
              {"sjet_2_pt", (nConstits > 1) ? tagged_constituents[1].pt()/MeV : 0. },

              //Third Jet
              {"sjet_3_mv2c10_binned", (nConstits > 2) ? get_mv2c10(tagged_constituents[2], jet_mv2c10_bands) : -1. },
              {"sjet_3_e", (nConstits > 2) ? tagged_constituents[2].E()/MeV : 0. },
              {"sjet_3_phi", (nConstits > 2) ? tagged_constituents[2].phi() : j.phi() },
              {"sjet_3_eta", (nConstits > 2) ? tagged_constituents[2].eta() : j.eta() },
              {"sjet_3_pt", (nConstits > 2) ? tagged_constituents[2].pt()/MeV : 0. }
            };

            outputs = _nn->compute(NN_Input);
            tag =  getTag(outputs);
            JetTags.push_back(tag);
          }

          const string category = get_3l_MCBOT_category(JetTags);
          _sigBins["3l_SR_"+category]->fill(); 
          _h["Hist3l_SR_"+category]->fill(HTjetlep);
        }

      }

      

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      
      //In validation mode, normalise histos to unity.
      if(_mode > 1){
        if (_h["DNN_V"]->integral() > 0){
          _h["DNN_V"]->normalize(1.);
        }
        if (_h["DNN_H"]->integral() > 0){
          _h["DNN_H"]->normalize(1.);
        }
        if (_h["DNN_top"]->integral() > 0){
          _h["DNN_top"]->normalize(1.);
        }
        if (_mode == 3){
          if (_h["DNN_light"]->integral() > 0){
            _h["DNN_light"]->normalize(1.);
          }
        }
        else {
          for (const string &s : {"DNN_V_forH", "DNN_V_fortop", "DNN_H_forV", "DNN_H_fortop",
                                  "DNN_top_forV", "DNN_top_forH", "DNN_light_forV", "DNN_light_forH",
                                  "DNN_light_fortop"}){
            if (_h[s]->integral() > 0){
            _h[s]->normalize(1.);
            }
          }
        }
        return;
      }

      const double sf = 1000*crossSection()*(139)/sumW();

      for (const string &s : _histoNames){
        _h[s]->scaleW(sf);
      }
      for (const string &s : _sigRegionNames){
        _sigBins[s]->scaleW(sf);
      }
      for (const string &s : _controlRegionNames){
        _controlBins[s]->scaleW(sf);
      }

      for (const string& s : {"2lEtmiss", "2lHTjet", "2lHTjetcentral", "2lEtmiss_HTjet", "2lEtmiss_HTjetcentral"}){
        _h[s]->scaleW(sf);
      }
    }

    /// @}


    private:
    /// @name Utility functions

    //trims reclustered jets by removing all subjets with pt < minpt. Does not work in place,
    //has seperate inout&output. Also returns a vector of new constituents so we definitely save this
    // info before the cluster sequence pointer disappears in a puff of smoke.
    //If I can use a fastjet trimmer to do this more elegantly please let me know.
    static void naive_trimmer(const Rivet::PseudoJets& input, Rivet::PseudoJets& output,
                       double pt_frac, vector<Rivet::PseudoJets>& new_constituents){
      output.clear();
      for (const PseudoJet& j : input){
        std::vector<size_t> to_keep;
        double min_pt = j.pt() * pt_frac;
        for (size_t i = 0; i < j.constituents().size(); ++i){
          if (j.constituents()[i].pt() >= min_pt){
            to_keep.push_back(i);
          }
        }
        //TODO this is overcomplicated.
        
        if (to_keep.size() > 0){
          PseudoJet jet;
          PseudoJets constits;
          for (const size_t i : to_keep){
            jet+=j.constituents()[i];
            constits.push_back(j.constituents()[i]);
          }
          output.push_back(jet);
          new_constituents.push_back(constits);
        }
      }
    }

    MCBot_TagType getTag(const map<string, double> &outputs){
      
      const bool hasVtag = outputs.at("dnnOutput_V") > _threshold.at("Vector");
      const bool hasHtag = outputs.at("dnnOutput_H") > _threshold.at("Higgs");
      const bool hasttag = outputs.at("dnnOutput_top") > _threshold.at("top");

      //No tags => Backgroundjet
      if (!(hasVtag || hasHtag || hasttag))
        return MCBot_TagType::Bkg;

      //Each of the single tag cases
      if (hasVtag && !(hasHtag || hasttag))
        return MCBot_TagType::Vec;
      if (hasHtag && !(hasVtag || hasttag))
        return MCBot_TagType::Higgs;
      if (hasttag && !(hasVtag || hasHtag))
        return MCBot_TagType::top;

      //Cases with two or three tags - this requires more work.
      // Pick the cas with highest individual score.
      // TODO: Technically speaking there's a pathlogical case where top-score is highest but it isn't top-tagged
      if (outputs.at("dnnOutput_V") > outputs.at("dnnOutput_H") &&
            outputs.at("dnnOutput_V") > outputs.at("dnnOutput_top")){
              return MCBot_TagType::Vec;
      }
      if (outputs.at("dnnOutput_H") > outputs.at("dnnOutput_V") &&
            outputs.at("dnnOutput_H") > outputs.at("dnnOutput_top")){
              return MCBot_TagType::Higgs;
      }
      if (outputs.at("dnnOutput_top") > outputs.at("dnnOutput_H") &&
            outputs.at("dnnOutput_top") > outputs.at("dnnOutput_V")){
              return MCBot_TagType::top;
      }

      // ELSE - SHOULDN'T EVER REACH THIS (but better safe than sorry
      // and it shuts the compiler up):
      MSG_ERROR("MCBot tagging not succesful -- debugging required");
      throw Error("MCBot tagging not succesful -- debugging required");
    }

    //Attempt at consistent btagging across working points.
    //Returns 100 if a jet isn't ever btagged,
    //or else the mv2c10 working point (60, 70, 77, 85)
    // at which the jet would be tagged,
    // TODO: Implement properly.
    int dummy_mv2c10_band(Jet& j){
      const double variate =  rand01();
      //TODO: mv2c10 cuts from ATL-PHYS-PUB-2016-012
      const double beff60 = JET_BTAG_EFFS(0.6, 1/34., 1/184., 1/1538.)(j);
      const double beff70 = JET_BTAG_EFFS(0.7, 1/12., 1/55., 1/381.)(j);
      const double beff77 = JET_BTAG_EFFS(0.77, 1/6., 1/22., 1/134.)(j);
      const double beff85 = JET_BTAG_EFFS(0.85, 1/3.1, 1/8.2, 1/33.)(j);
      //For the purposes of this analysis outside the NN, a b-tag is done at 77%
      const bool btag  = rand01() < beff77 && j.abseta() < 2.5;
      // Remove b-tags if needed, and add a dummy one if needed
      if (!btag && j.bTagged()) j.tags().erase(std::remove_if(j.tags().begin(), j.tags().end(), hasBottom), j.tags().end());
      if (btag && !j.bTagged()) j.tags().push_back(Particle(PID::BQUARK, j.mom()));
      
      //Ignore c-tagging.
      if (j.abseta() >= 2.5) return 100;
      else if (variate < beff60) return 60;
      else if (variate < beff70) return 70;
      else if (variate < beff77) return 77;
      else if (variate < beff85) return 85;
      else return 100;
    }

    static double mv2c10_score_from_band(const int band){
      switch (band){
        case(60):
          return 0.934906;
        case(70):
          return 0.8244273;
        case(77):
          return 0.645925;
        case(85):
          return 0.1758475;
        case(100):
          return 0.0;
        default:
          throw Error("Failed to convert mv2c10 band to score");
      }
    }

    // Utility function summing up previous two in one go.
    double get_mv2c10(const Jet &j, const vector<pair<ThreeMomentum, int>>& mv2c10_bands){
      auto it = std::find_if(mv2c10_bands.begin(), mv2c10_bands.end(),
                                  [&j](const pair<ThreeMomentum, int>& p){return (deltaR(p.first, j.p3()) < 0.1);});
      if (it != mv2c10_bands.end()){
        return mv2c10_score_from_band(it->second);
      }
      else {
        throw Error("FAILED TO FIND mv2c10");
      }
    }

    //Get the MCBot category (defined in Table 2) for 2l channel
    string get_2l_MCBOT_category(vector<MCBot_TagType> tags, size_t nbtags){
      size_t nVtags = 0, nHtags = 0, ntoptags = 0;
      for(const MCBot_TagType tt : tags){
        if (tt == MCBot_TagType::Vec)
          ++nVtags;
        else if (tt == MCBot_TagType::Higgs){
          ++nHtags;
        }
        else if (tt == MCBot_TagType::top){
          ++ntoptags;
        }
      }
      if (nVtags == 0 && nHtags == 0 && ntoptags == 0){
        return "Notag";
      }
      else if (nVtags == 1 && nHtags == 0 && ntoptags == 0){
        return "Vtag";
      }
      else if (nVtags == 0 && nHtags == 1 && ntoptags == 0){
        return "Htag";
      }
      else if (nVtags == 0 && nHtags == 0 && ntoptags == 1){
        return "toptag";
      }
      else if ((nVtags == 2 && nHtags == 0 && ntoptags == 0) ||
        (nVtags == 0 && nHtags == 2 && ntoptags == 0)
        || (nVtags == 1 && nHtags == 0 && ntoptags == 1 && nbtags == 1)
        || (nVtags == 1 && nHtags == 1 && ntoptags == 0 && nbtags >= 2)
        || (nVtags == 0 && nHtags == 0 && ntoptags == 2 && nbtags >= 2)){
        return "Doubletag1";
      }
      else if ((nVtags == 0 && nHtags == 1 && ntoptags == 1)
        || (nVtags == 0 && nHtags == 0 && ntoptags == 2 && nbtags == 1)){
        return "Doubletag2";
      }
      else if ((nVtags == 1 && nHtags == 1 && ntoptags == 0 && nbtags == 1)
        || (nVtags == 1 && nHtags == 0 && ntoptags == 1 && nbtags >= 2)
        || (nVtags + nHtags + ntoptags > 2)){
        return "OF";
      }
      else {
        MSG_ERROR("FAILED TO CATEGORISE JET - DEBUGGING REQUIRED");
        throw Error("FAILED TO CATEGORISE JET - DEBUGGING REQUIRED");
      }
    }


    //Get the MCBot category (defined in Table 2) for 2l channel
    string get_3l_MCBOT_category(vector<MCBot_TagType> tags){
      size_t nVtags = 0, nHtags = 0, ntoptags = 0;
      for(const MCBot_TagType tt : tags){
        if (tt == MCBot_TagType::Vec)
          ++nVtags;
        else if (tt == MCBot_TagType::Higgs){
          ++nHtags;
        }
        else if (tt == MCBot_TagType::top){
          ++ntoptags;
        }
      }
      if (nVtags == 0 && nHtags == 0 && ntoptags == 0){
        return "Notag";
      }
      else if (nVtags >= 1 && nHtags == 0 && ntoptags == 0){
        return "Vtag";
      }
      else if (nVtags == 0 && nHtags >= 1 && ntoptags == 0){
        return "Htag";
      }
      else if (nVtags == 0 && nHtags == 0 && ntoptags >= 1){
        return "toptag";
      }
      else 
        return "OF";
    }

    //Very rough and hacky
    //Find all Vector bosons, higgs, and tops that are not children of themselves.
    //TODO: implement as a projection?
    static void getVHandtop_fromEvent(std::vector<Particle>& ps, const Event& ge, 
                                      std::vector<int> wanted_pids = {6,23,24,25}){
      for (const auto& p : ge.allParticles()){
        if (p.pt() < 20*GeV){
          continue;
        }
        const int pid = p.pid();
        //TODO: Should there also be a status code check?
        //If its Vht, check its parents aren't the same particle to avoid including what is really the same particle many times.
        for (const int targetPID : wanted_pids){
          if (abs(pid) == targetPID){
            Particles parents = p.parents();
            auto it = std::find_if(parents.begin(), parents.end(),
                                  [pid](const Particle &p){return p.pid() == pid;});
            if (it == parents.end()){
              ps.push_back(p);
            }
          }
        }
      }
    }

    /// }

    /// @name Member variables

    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, CounterPtr> _sigBins;
    map<string, CounterPtr> _controlBins;
    /// @}

    int _mode;
    std::unique_ptr<lwt::LightweightNeuralNetwork> _nn;
    const std::map<string, double> _threshold = {{"Vector", 0.26}, {"Higgs", 0.31}, {"top", 0.21}};

    //Some names for easy iteration and cleaner, shorter code
    const vector<string> _2L_MCBOT_categories = {"Notag", "Vtag", "Htag", "toptag", "Doubletag1", "Doubletag2", "OF"};
    const vector<string> _3L_MCBOT_categories = {"Notag", "Vtag", "Htag", "toptag", "OF"};
    const vector<string> _histoNames = {"Hist2l_1b_CR", "Hist2l_2b_CR",
                                        "Hist2l_2b_SR_OF", "Hist2l_2b_SR_Doubletag1", "Hist2l_2b_SR_Doubletag2",
                                        "Hist2l_2b_SR_Vtag", "Hist2l_2b_SR_Htag", "Hist2l_2b_SR_toptag", "Hist2l_2b_SR_Notag",
                                        "Hist2l_1b_SR_OF", "Hist2l_1b_SR_Doubletag1", "Hist2l_1b_SR_Doubletag2",
                                        "Hist2l_1b_SR_Vtag", "Hist2l_1b_SR_Htag", "Hist2l_1b_SR_toptag", "Hist2l_1b_SR_Notag",
                                        "Hist3l_VV_CR", "Hist3l_SR_Notag", "Hist3l_SR_Vtag", "Hist3l_SR_Htag",
                                        "Hist3l_SR_toptag", "Hist3l_SR_OF"};
    const vector<string> _sigRegionNames = {"2l_2b_SR_OF", "2l_2b_SR_Doubletag1", "2l_2b_SR_Doubletag2",
                                        "2l_2b_SR_Vtag", "2l_2b_SR_Htag", "2l_2b_SR_toptag", "2l_2b_SR_Notag",
                                        "2l_1b_SR_OF", "2l_1b_SR_Doubletag1", "2l_1b_SR_Doubletag2",
                                        "2l_1b_SR_Vtag", "2l_1b_SR_Htag", "2l_1b_SR_toptag", "2l_1b_SR_Notag",
                                        "3l_SR_Notag", "3l_SR_Vtag", "3l_SR_Htag", "3l_SR_toptag", "3l_SR_OF"};
    const vector<string> _controlRegionNames = {"2l_1b_CR", "2l_2b_CR","3l_VV_CR"};




    /// @}


  };


  RIVET_DECLARE_PLUGIN(ATLAS_2022_I2172216);

}