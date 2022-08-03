// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/SmearedJets.hh"
#include "Rivet/Projections/SmearedParticles.hh"
#include "Rivet/Projections/SmearedMET.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Tools/MCBot_tagger.hh"
#include "Rivet/Tools/RivetPaths.hh"
#include "Rivet/Tools/Random.hh"
#include "Rivet/Tools/JetSmearingFunctions.hh"
#include "Rivet/Math/MathUtils.hh"
#include "fastjet/contrib/VariableRPlugin.hh"
#include "fastjet/tools/Filter.hh"
#include <fstream>

namespace Rivet {

//trims reclustered jets by removing all subjets with pt < minpt. Does not work in place,
//has seperate inout&output. Also returns a vector of new constituents so we definitely save this
// info before the cluster sequence pointer disappears in a puff of smoke.
//If I can use a fastjet trimmer to do this more elegantly please let me know.
void naive_RC_trimmer(const Rivet::PseudoJets& input, Rivet::PseudoJets& output,
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

//Numbers match MCBot_tagtype.
//TODO: probably should only have one type of enum
enum class DNN_Category{
  V_tag=2,
  H_tag=3,
  top_tag=1,
  bkg_tag=0
};

    /// @brief Add a short analysis description here
  class ATLAS_2018_I1685207 : public Analysis {


  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2018_I1685207);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // BCG = background mode, SGN = signal mode (V, H, top)  
      // default to signal mode

      // electrons
      const PromptFinalState fse(Cuts::abseta < 2.47 && Cuts::pt > 20*GeV && 
		                             Cuts::abspid == PID::ELECTRON &&
								                 !Cuts::absetaIn(1.37, 1.52));
	    declare("Elec", fse);
	    SmearedParticles recoelectrons(fse, ELECTRON_RECOEFF_ATLAS_RUN2, ELECTRON_SMEAR_ATLAS_RUN2);
    	declare(recoelectrons, "SmearedElec");

	    /// muons
	    const PromptFinalState fsm(Cuts::abseta < 2.5 && Cuts::pT > 20*GeV && Cuts::abspid == PID::MUON);
	    declare("Muon", fsm);
      //TODO: Check this is the right smearing.
	    SmearedParticles recomuons(fsm, MUON_EFF_ATLAS_RUN2, MUON_SMEAR_ATLAS_RUN2);
    	declare(recomuons, "SmearedMuon");

	    /// jets
	    const FinalState fsj(Cuts::abseta < 4.8);
	    
	    FastJets Sj(fsj, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE);
	    declare("Sjet", Sj);
      
      // No smearing:
	    //SmearedJets SSj(Sj, JET_SMEAR_IDENTITY, JET_BTAG_EFFS(0.77, 1./6.2, 1./134));
      /// Energy-resolution smearing only:
	    SmearedJets SSj(Sj, JET_SMEAR_ATLAS_RUN2, JET_BTAG_EFFS(0.77, 1./6.2, 1./134));
      /// @todo Also look into angular smearing? Need a custom smearing function, building on the ATLAS R2
      //SmearedJets SSj(Sj, JET_SMEAR_ANGULAR, JET_BTAG_EFFS(0.77, 1./6.2, 1./134));
    	declare(SSj, "smearedSjet");
      
      // pTmiss
      declare(VisibleFinalState(Cuts::abseta < 4.8),"vfs");

	    // Book histograms and counters
      //n.b. the bins values, expected etc are given in the paper, but are NOT in HEPDATA.
      //Signal bins
      book(_sigBins["VV_0t_2b"], "VV_0t_2b");
      book(_sigBins["VV_0t_3b"], "VV_0t_3b");
      book(_sigBins["VV_1t_2b"], "VV_1t_2b");
      book(_sigBins["VV_1t_3b"], "VV_1t_3b");

      book(_sigBins["VH_0t_2b"], "VH_0t_2b");
      book(_sigBins["VH_0t_3b"], "VH_0t_3b");
      book(_sigBins["VH_1t_2b"], "VH_1t_2b");
      book(_sigBins["VH_1t_3b"], "VH_1t_3b");

      book(_sigBins["HH_0t_3b"], "HH_0t_3b");
      book(_sigBins["HH_1t_3b"], "HH_1t_3b");
      book(_sigBins["XX_2t_2b"], "XX_2t_2b");
      book(_sigBins["XX_2t_3b"], "XX_2t_3b");

      //Validation region bins
      book(_valBins["HH_0t_2b"], "HH_0t_2b");
      book(_valBins["HH_1t_2b"], "HH_1t_2b");

      book(_valBins["VV_0t_1b"], "VV_0t_1b");
      book(_valBins["VV_1t_1b"], "VV_1t_1b");
      book(_valBins["VH_0t_1b"], "VH_0t_1b");
      book(_valBins["VH_1t_1b"], "VH_1t_1b");
      book(_valBins["HH_0t_1b"], "HH_0t_1b");
      book(_valBins["HH_1t_1b"], "HH_1t_1b");
      book(_valBins["XX_2t_1b"], "XX_2t_1b");

      //Find the json file
      //const std::string nn_datafilename = "ATLAS_2018_I1685207.nn.json.yoda";
      const std::string nn_datafilename = "ATLAS_2018_I1685207.btagEnriched.nn.json.yoda";
      //TODO: Would be nice to use the proper find syntax but there seems to be assumptions
      // about .yoda endings. Someone who understands the paths system better would do it more
      // elegantly.
      const std::string nn_datafilepath =  getDataPath()+"/Rivet/"+nn_datafilename;
      
      _MCbottagger = std::make_unique<MCBot_tagger>(MCBot_tagger(nn_datafilepath));

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      Particles smeared_electrons = apply<ParticleFinder>(event, "SmearedElec").particles();
      //This is a 0-lepton analysis.
      if (smeared_electrons.size() != 0){
        vetoEvent;
      }
      Particles smeared_muons = apply<ParticleFinder>(event, "SmearedMuon").particles();
      //This is a 0-lepton analysis.
      if (smeared_muons.size() != 0){
        vetoEvent;
      }
      
      Jets smeared_small_jets = apply<JetAlg>(event, "smearedSjet").jetsByPt();
      idiscard(smeared_small_jets, Cuts::pt <= 25*GeV && Cuts::abseta >= 2.5); 
      
      //Todo JVT means we lose 8% of small jets?
      //Get b-tagged jets:
      const Jets smeared_bjets = filter_select(smeared_small_jets, hasBTag());

      FourMomentum pTmiss;
      for (const Particle& p : apply<VisibleFinalState>(event, "vfs").particles() ) {
        pTmiss -= p.momentum();
      }
      double ETmiss = pTmiss.pT();
      

      //Recluster small jets using variable R
      // rho= 315GeV, 0.4<Reff<1.2
      //TODO: probably more efficient to move after pre-selection.
      fastjet::JetDefinition::Plugin* variableRplugin = new fastjet::contrib::VariableRPlugin(315*GeV, 0.4, 1.2, fastjet::contrib::VariableRPlugin::AKTLIKE);
      fastjet::JetDefinition jdef(variableRplugin);
      ClusterSequence cseq(smeared_small_jets, jdef);
      PseudoJets VRC_jets = cseq.inclusive_jets();
    
      // Trimming
      PseudoJets TrimmedVRCjets;
      vector<PseudoJets> NewConstits;
      naive_RC_trimmer(VRC_jets, TrimmedVRCjets, 0.05, NewConstits);
      
      // pT, eta, and mass cuts
      PseudoJets FilteredVRCjets;
      vector<PseudoJets> FilteredNewConstits;

      PseudoJet tr_pj;
      for (size_t i = 0; i < TrimmedVRCjets.size(); ++i){
        tr_pj = TrimmedVRCjets[i];
        if (tr_pj.pt() > 150*GeV && abs(tr_pj.eta()) < 2.5 && tr_pj.m() > 40*GeV){
          FilteredVRCjets.push_back(tr_pj);
          FilteredNewConstits.push_back(NewConstits[i]);
        }
      }
      
      int VRCsize = FilteredVRCjets.size();
      if (VRCsize == 0) {
        vetoEvent;
      }
      
      // signal - vRC jets for further analysis with their corresponding constituent jets in SignalConstits
      PseudoJets signal; 
      vector<PseudoJets> SignalConstits;
      
      signal = FilteredVRCjets;
      SignalConstits = FilteredNewConstits;

      // DNN tags of vRC jets
      std::vector<DNN_Category> VRCjet_tags;
      
      // DNN scores (probabilities)
      std::map<string, double> outputs;

      // analysis of selected vRC-jets
      for (size_t counter = 0; counter < signal.size(); ++counter){
        auto& j = signal[counter];
          
        //First we need a Jets of all the constituents, that has b-tagging info. 
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
        
        //computes the D values (probabilities)
        _MCbottagger->computeScores(j, tagged_constituents, outputs);

        //Tag the jets in an approximation of the MCBot NN.
        MCBot_TagType DNNtag = _MCbottagger->tag(j,tagged_constituents);
        VRCjet_tags.push_back(static_cast<DNN_Category>(DNNtag));
      }
      
      //Preselection
      //HT > 1250GeV
      // HT = total scalar sum of the transverse momenta of all track particles and energy deposits
      double HT = std::accumulate(smeared_small_jets.begin(), smeared_small_jets.end(), 0*GeV, 
                                  [](double count, const Jet& nextjet){return (count+nextjet.pt());}); 
      if (HT <= 1250*GeV){
        vetoEvent;
      }

      //ETMiss < 200GeV and ETmiss > 40GeV
      if (ETmiss <= 40*GeV || ETmiss >= 200*GeV){
        vetoEvent;
      }
      
      //Four small jets 
      if (smeared_small_jets.size() < 4){
        vetoEvent;
      }
      
      //Minimum pT's for leading smalljets 
      if (smeared_small_jets[0].pT() <= 300*GeV || smeared_small_jets[1].pT() <= 200*GeV ||
        smeared_small_jets[2].pT() <= 125*GeV || smeared_small_jets[3].pT() <= 75*GeV){
        vetoEvent;
      }
    
      //At least 2 bjets
      if (smeared_bjets.size() < 2){
        vetoEvent;
      }
      
      //Two vRC jets tagged V or H
      int nVtags = std::count_if(VRCjet_tags.begin(), VRCjet_tags.end(), 
                                  [](const DNN_Category dc){return (dc == DNN_Category::V_tag);});
      int nHtags = std::count_if(VRCjet_tags.begin(), VRCjet_tags.end(), 
                                  [](const DNN_Category dc){return (dc == DNN_Category::H_tag);});
      int ntoptags = std::count_if(VRCjet_tags.begin(), VRCjet_tags.end(), 
                                  [](const DNN_Category dc){return (dc == DNN_Category::top_tag);});                              

      if (nVtags + nHtags < 2){
        vetoEvent;
      }
      
      //select signal/validation/control region.
      //TODO: The pre-selection cuts say 2 or more (v or H) tagged jets, but each signal region requires only two.
      // I also can't see some sort of tie-break procedure outlined (e.g. take tags of two highest pT jets)
      // Unless this is what the XX bins are for?
      //Only one of SR and CR is filled (please!)
      //VV signal regions
      if (nVtags == 2 && nHtags == 0 && ntoptags == 0 && smeared_bjets.size() == 2){
        _sigBins["VV_0t_2b"]->fill();
      }
      else if (nVtags == 2 && nHtags == 0 && ntoptags == 0 && smeared_bjets.size() > 2){
        _sigBins["VV_0t_3b"]->fill();
      }
      else if (nVtags == 2 && nHtags == 0 && ntoptags == 1 && smeared_bjets.size() == 2) {
        _sigBins["VV_1t_2b"]->fill();
      } 
      else if (nVtags == 2 && nHtags == 0 && ntoptags == 1 && smeared_bjets.size() > 2) {
        _sigBins["VV_1t_3b"]->fill();
      }
      //VH signal regions
      else if (nVtags == 1 && nHtags == 1 && ntoptags == 0  && smeared_bjets.size() == 2){
        _sigBins["VH_0t_2b"]->fill();
      }
      else if (nVtags == 1 && nHtags == 1 && ntoptags == 0 && smeared_bjets.size() > 2) {
        _sigBins["VH_0t_3b"]->fill();
      }
      else if (nVtags == 1 && nHtags == 1 && ntoptags == 1 && smeared_bjets.size() == 2){
        _sigBins["VH_1t_2b"]->fill();
      }
      else if (nVtags == 1 && nHtags == 1 && ntoptags == 1 && smeared_bjets.size() > 2 ) {
        _sigBins["VH_1t_3b"]->fill();
      }
      //HH signal regions
      else if (nVtags == 0 && nHtags == 2 && ntoptags == 0 && smeared_bjets.size() > 2){
        _sigBins["HH_0t_3b"]->fill();
      }
      else if ( nVtags == 0 && nHtags == 2 && ntoptags == 1 && smeared_bjets.size() > 2){
        _sigBins["HH_1t_3b"]->fill();
      }
      //XX signal regions.
      //TODO: mild inconcistency between text and table: try both == and >= and see which match better the results...
      
      else if (nVtags + nHtags >= 2 && ntoptags >= 2 && smeared_bjets.size()==2){
        _sigBins["XX_2t_2b"]->fill();
      }
      else if (nVtags + nHtags >= 2 && ntoptags >= 2 && smeared_bjets.size()>=3){
        _sigBins["XX_2t_3b"]->fill();
      }

      //Now cover the validation regions

      //MultiJet Background Uncertainty Validation Regions (2)
      //TODO: This makes no sense. They say there are two regions but define one?!?
      //Possibly involves altering nVtags? But how?!
      //If I had to guess, by symmetry its HH_1t_2b, but this contradicts the paper. Going with it for now...
      else if (smeared_bjets.size()==2 && nHtags == 2 && ntoptags == 0 && nVtags == 0){
        _valBins["HH_0t_2b"]->fill();
      }
      else if (smeared_bjets.size()==2 && nHtags == 2 && ntoptags == 1 && nVtags == 0){
        _valBins["HH_1t_2b"]->fill();
      }
      //Closure uncertainty validation regions (7)
      else if (nVtags == 2 && nHtags == 0 && ntoptags == 0 && smeared_bjets.size() == 1){
        _valBins["VV_0t_1b"]->fill();
      }
      else if (nVtags == 2 && nHtags == 0 && ntoptags == 1 && smeared_bjets.size() == 1) {
        _valBins["VV_1t_1b"]->fill();
      } 
      else if (nVtags == 1 && nHtags == 1 && ntoptags == 0  && smeared_bjets.size() == 1){
        _valBins["VH_0t_1b"]->fill();
      } 
      else if (nVtags == 1 && nHtags == 1 && ntoptags == 1 && smeared_bjets.size() == 1){
        _valBins["VH_1t_1b"]->fill();
      } 
      else if (nVtags == 0 && nHtags == 2 && ntoptags == 0 && smeared_bjets.size() == 1){
        _valBins["HH_0t_1b"]->fill();
      } 
      else if ( nVtags == 0 && nHtags == 2 && ntoptags == 1 && smeared_bjets.size() == 1){
        _valBins["HH_1t_1b"]->fill();
      }
      else if (nVtags + nHtags >= 2 && ntoptags >= 2 && smeared_bjets.size()==1){
        _valBins["XX_2t_1b"]->fill();
      }

      //DEBUG ONLY
      // else {
      //   std::cout << "WARNING, event passed preselection but belonged in no bins\n";
      //   std::cout << "V, H, t, b = (" << nVtags << ", " << nHtags << ", " << ntoptags << ", " << smeared_bjets.size() << ")\n";
      // }
    }
      

      /// Normalise histograms etc., after the run
      void finalize() {

        //Shortcut for validation only:
        std::cout << "\nSIG BINS\n";
        std::cout << _sigBins["VV_0t_2b"]->val() << ", " << _sigBins["VV_0t_3b"]->val() << ", ";
        std::cout << _sigBins["VV_1t_2b"]->val() << ", " << _sigBins["VV_1t_3b"]->val() << ", ";
        std::cout << _sigBins["VH_0t_2b"]->val() << ", " << _sigBins["VH_0t_3b"]->val() << ", ";
        std::cout << _sigBins["VH_1t_2b"]->val() << ", " << _sigBins["VH_1t_3b"]->val() << ", ";
        std::cout << _sigBins["HH_0t_3b"]->val() << ", " << _sigBins["HH_1t_3b"]->val() << ", ";
        std::cout << _sigBins["XX_2t_2b"]->val() << ", " << _sigBins["XX_2t_3b"]->val() << ", ";

        // scale(_h, crossSectionPerEvent() / femtobarn);
        // std::cout << "\nSIG BINS (scaled to XS)\n";
        // std::cout << _sigBins["VV_0t_2b"]->val() << ", " << _sigBins["VV_0t_3b"]->val() << ", ";
        // std::cout << _sigBins["VV_1t_2b"]->val() << ", " << _sigBins["VV_1t_3b"]->val() << ", ";
        // std::cout << _sigBins["VH_0t_2b"]->val() << ", " << _sigBins["VH_0t_3b"]->val() << ", ";
        // std::cout << _sigBins["VH_1t_2b"]->val() << ", " << _sigBins["VH_1t_3b"]->val() << ", ";
        // std::cout << _sigBins["HH_0t_3b"]->val() << ", " << _sigBins["HH_1t_3b"]->val() << ", ";
        // std::cout << _sigBins["XX_2t_2b"]->val() << ", " << _sigBins["XX_2t_3b"]->val() << ", ";
      }

      /// @}


      /// @name Histograms and Counters
      /// @{
      map<string, CounterPtr> _sigBins;
      map<string, CounterPtr> _valBins;
      map<string, Histo1DPtr> _h;
      
      
      std::unique_ptr<MCBot_tagger> _MCbottagger;

    };


  RIVET_DECLARE_PLUGIN(ATLAS_2018_I1685207);


}
