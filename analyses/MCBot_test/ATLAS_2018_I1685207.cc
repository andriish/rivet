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
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/tools/Filter.hh"
#include <fstream>

//trims reclustered jets by removing all subjets with pt < minpt. Does not work in place,
//has seperate inout&output. Also returns a vector of new constituents so we definitely save this
// info before the cluster sequence pointer disappears in a puff of smoke.
//If I can use a fastjet trimmer to do this more elegantly please let me know.
void naive_RC_trimmer(const Rivet::PseudoJets& input, Rivet::PseudoJets& output,
                       double pt_frac, vector<Rivet::PseudoJets>& new_constituents){
  output.clear();
  for (const Rivet::PseudoJet& j : input){
    std::vector<size_t> to_keep;
    double min_pt = j.pt() * pt_frac;
    for (size_t i = 0; i < j.constituents().size(); ++i){
      if (j.constituents()[i].pt() >= min_pt){
        to_keep.push_back(i);
      }
    }
    //TODO this is overcomplicated.
    
    if (to_keep.size() > 0){
      Rivet::PseudoJet jet;
      Rivet::PseudoJets constits;
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
//TOOD: probably should only have one type of enum
enum class DNN_Category{
  V_tag=2,
  H_tag=3,
  top_tag=1,
  bkg_tag=0
};

//Very rough and hacky
//Find all Vector bosons, higgs, and tops that are not children of themselves.
//TODO: implement as a projection?
void getVHandtop_fromEvent(std::vector<Rivet::Particle>& ps, const Rivet::Event& ge){
  for (const auto& p : ge.allParticles()){
    int pid = p.pid();
    //TODO: Should there also be a status code check?
    //If its Vht, check its parents aren't the same particle to avoid including what is really the same particle many times.
    if (abs(pid) == 23 || abs(pid) == 24 || abs(pid) == 6 || abs(pid) == 25){
      Rivet::Particles parents = p.parents();
      auto it = std::find_if(parents.begin(), parents.end(),
                             [pid](const Rivet::Particle &p){return p.pid() == pid;});
      if (it == parents.end()){
        ps.push_back(p);
      }
    }
  }

}

// 
DNN_Category getTrueDNNtag(const Rivet::PseudoJet& pj, const Rivet::Particles& comp){
  for (const Rivet::Particle& comp_p : comp){
    //Jets need to be within deltaR 0.75rho/pT
    //TODO: I guess rho is the same from vRC jet reclustering but is it the particle or jet pT?
    // I think jet but not 100% sure.
    //TODO: In the paper when they ran the simulation to train the DNN, they discarded events with
    // with multiple "true" matches. But what is correct in this context.
    // For now I'm going with quick and dirty "first one in range works";
    if (Rivet::deltaR(Rivet::momentum3(pj), Rivet::momentum3(comp_p.pseudojet())) < 315*Rivet::GeV/pj.pt()){
      if (comp_p.abspid() == 23 || comp_p.abspid() == 24){
        //std::cout << "TAGGING Vector Boson" << std::endl;
        return DNN_Category::V_tag;
      }
      if (comp_p.abspid() == 25){
        ///std::cout << "TAGGING Higgs" << std::endl;
        return DNN_Category::H_tag;
      }
      if (comp_p.abspid() == 6){
        ///std::cout << "TAGGING top" << std::endl;
        return DNN_Category::top_tag;
      }
    }
  }

  //std::cout << "TAGGING bkg" << std::endl;
  return DNN_Category::bkg_tag;
}


// angular smearing function: building on JET_SMEAR_ATLAS_RUN2 (copy of run 1)
namespace Rivet {

Jet JET_SMEAR_ANGULAR2(const Jet& j) {
  // Jet energy resolution lookup
  //   original -- Implemented by Matthias Danninger for GAMBIT, based roughly on
  //   https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2015-017/
  //   Parameterisation can be still improved, but eta dependence is minimal
  /// @todo Also need a JES uncertainty component?
  static const vector<double> binedges_pt = {0., 50., 70., 100., 150., 200., 1000., 10000.};
  static const vector<double> jer = {0.145, 0.115, 0.095, 0.075, 0.07, 0.05, 0.04, 0.04}; //< note overflow value
  const int ipt = binIndex(j.pt()/GeV, binedges_pt, true);
  if (ipt < 0) return j;
  const double resolution = jer.at(ipt);

  // Smear by a Gaussian centered on 1 with width given by the (fractional) resolution
  /// @todo Is this the best way to smear? Should we preserve the energy, or pT, or direction?
  const double fsmear = max(randnorm(1., resolution), 0.); // use this for smearing the energy...
  const double mass = j.mass2() > 0 ? j.mass() : 0; //< numerical carefulness... is it correct to calculate the mass here? mass vs smeared energy inconsistency???
  
  Jet j1(Rivet::FourMomentum::mkXYZM(j.px()*fsmear, j.py()*fsmear, j.pz()*fsmear, mass));

  // smearing in eta-phi ...

  double dsmear = max(randnorm(0., 0.1), 0.);
  double theta = rand01() * M_2_PI; // why 2/pi instead of 2pi??
  // double theta = 2* rand01() * M_PI; // alternatively??
    
  
  return Jet(FourMomentum::mkEtaPhiME(j.eta()+dsmear*cos(theta), mapAngle0To2Pi(j.phi()+dsmear*sin(theta)), j1.mass(), j1.E()));  
}

Jet JET_SMEAR_ANGULAR3(const Jet& j) {
  // Jet energy resolution lookup
  //   original -- Implemented by Matthias Danninger for GAMBIT, based roughly on
  //   https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2015-017/
  //   Parameterisation can be still improved, but eta dependence is minimal
  /// @todo Also need a JES uncertainty component?
  static const vector<double> binedges_pt = {0., 50., 70., 100., 150., 200., 1000., 10000.};
  static const vector<double> jer = {0.145, 0.115, 0.095, 0.075, 0.07, 0.05, 0.04, 0.04}; //< note overflow value
  const int ipt = binIndex(j.pt()/GeV, binedges_pt, true);
  if (ipt < 0) return j;
  const double resolution = jer.at(ipt);

  // Smear by a Gaussian centered on 1 with width given by the (fractional) resolution
  /// @todo Is this the best way to smear? Should we preserve the energy, or pT, or direction?
  const double fsmear = max(randnorm(1., resolution), 0.); // use this for smearing the energy...
  const double mass = j.mass2() > 0 ? j.mass() : 0; //< numerical carefulness... is it correct to calculate the mass here? mass vs smeared energy inconsistency???
  
  Jet j1(Rivet::FourMomentum::mkXYZM(j.px()*fsmear, j.py()*fsmear, j.pz()*fsmear, mass));

  // smearing in eta-phi ...

  double dsmear = max(randnorm(0., 1.), 0.);
  double theta = rand01() * M_2_PI; // why 2/pi instead of 2pi??
  // double theta = 2* rand01() * M_PI; // alternatively??
    
  
  return Jet(FourMomentum::mkEtaPhiME(j.eta()+dsmear*cos(theta), mapAngle0To2Pi(j.phi()+dsmear*sin(theta)), j1.mass(), j1.E()));  
}


    /// @brief Add a short analysis description here
  class ATLAS_2018_I1685207 : public Analysis {


  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2018_I1685207);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // BCG = background mode, SGN = signal mode   
      // default to signal mode

      _mode = 1; 
      if (getOption("MODE") == "BCG"){
        _mode = 0;
      }
      else if (getOption("MODE") == "SGN") {
        _mode = 1;
      }

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
	    //SmearedJets SSj(Sj, JET_SMEAR_ATLAS_RUN2, JET_BTAG_EFFS(0.77, 1./6.2, 1./134));
      /// @todo Also look into angular smearing? Need a custom smearing function, building on the ATLAS R2
      SmearedJets SSj(Sj, JET_SMEAR_ANGULAR3, JET_BTAG_EFFS(0.77, 1./6.2, 1./134));
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

      // discriminant function histograms
      book(_h["PV"], "PV",45,-3,1.5);
      book(_h["PH"], "PH",55,-3,2.5);
      book(_h["Pt"], "Pt",55,-3,2.5);      

      // Les Houches angularity, number of constituent small-R jets
      book(_h["LHA"],"LHA",200, 0, 1);
      book(_h["Nconsts"],"Nconsts", 10,0,10);
      // book(_h["tau21"], "tau21", 100, 0, 1);
      // book(_h["tau32"], "tau32", 100, 0, 1);
      

      //Find the json file
      const std::string nn_datafilename = "ATLAS_2018_I1685207.nn.json.yoda";
      //TODO: Would be nice to use the proper find syntax but there seems to be assumptions
      // about .yoda endings. Someone who understands the paths system better would do it more
      // elegantly.
      const std::string nn_datafilepath =  getDataPath()+"/Rivet/"+nn_datafilename;
      
      _MCbottagger = std::make_unique<MCBot_tagger>(MCBot_tagger(nn_datafilepath));

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static int count = 0;
      count++;

      Particles smeared_electrons = apply<ParticleFinder>(event, "SmearedElec").particles();
      //This is a 0-lepton analysis.
      if (smeared_electrons.size() != 0){
        vetoEvent;
      }
      Particles smeared_muons = apply<ParticleFinder>(event, "SmearedMuon").particles();
      //This is a 0=lepton analysis.
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
      
      
      //TODO there's got to be an easier way than this mess.
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
      
      // defining signal as an input for further analysis and filling signal constits with corresponding constituent jets
      // for background mode: analyse all TrimmedVRCjets, for signal mode, pick the 2 highest pT and then analyse the one with higher mass

      PseudoJets signal; 
      vector<PseudoJets> SignalConstits;

      // get the V, H, and tops from the event
      std::vector<Particle> VHandtops;
      getVHandtop_fromEvent(VHandtops, event); 
      
      
      if (_mode == 0){

        // input variables for the for the loop
        double pTjet = 0;
        double largest = 0;
        int index = 0;
        
        for (int i = 0; i < VRCsize;++i){
          pTjet = FilteredVRCjets[i].pt();
          
          if (pTjet >= largest){
            largest = pTjet;
            index = i;
            }
        }

        signal.push_back(FilteredVRCjets[index]);
        SignalConstits.push_back(FilteredNewConstits[index]);
      }
      
      else if (_mode  == 1) {

        // input variables for the for the loop
        double pTjet = 0;
        double largest1 = 0;
        double largest2 = 0;
        int index1 = 0;
        int index2 = 0;

        for (int i = 0; i < VRCsize;++i){
          pTjet = FilteredVRCjets[i].pt();
          
          if (pTjet >= largest1){
            largest2 = largest1;
            largest1 = pTjet;
            index2 = index1;
            index1 = i;
            }
          else if (pTjet >= largest2){
            largest2 = pTjet;
            index2 = i;
          }

        }

        auto it_1 = std::find_if(VHandtops.begin(), VHandtops.end(),
                                  [&FilteredVRCjets,index1](const Particle& VHTop){return (deltaR(momentum3(VHTop.pseudojet()), momentum3(FilteredVRCjets[index1])) < 0.1);});

        auto it_2 = std::find_if(VHandtops.begin(), VHandtops.end(),
                                  [&FilteredVRCjets,index2](const Particle& VHTop){return (deltaR(momentum3(VHTop.pseudojet()), momentum3(FilteredVRCjets[index2])) < 0.1);});


        if (FilteredVRCjets[index1].m() >= FilteredVRCjets[index2].m() && it_1 != VHandtops.end()){ 
            signal.push_back(FilteredVRCjets[index1]);
            SignalConstits.push_back(FilteredNewConstits[index1]);
        }
        else if (it_2 != VHandtops.end()){
          signal.push_back(FilteredVRCjets[index2]);
          SignalConstits.push_back(FilteredNewConstits[index2]);
        }
        else {
          vetoEvent;
        }

      }

      

      
      // DNN tags of vRC jets
      std::vector<DNN_Category> VRCjet_tags;
      
      // DNN scores (probabilities) -- output of the NN 
      std::map<string, double> outputs;
      
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
        

        // Les Houches Angularity
        
        double LHA_sum = 0;
        
        for (const PseudoJet& pj : SignalConstits[counter]){
          double pT = pj.pt();
          double dR = pj.squared_distance(j);
          LHA_sum += (pT * pow(dR, 0.25)); 
        }
        
        double LHA = LHA_sum/j.pt();
        
        _h["LHA"]->fill(LHA);
        
        // number of constituent jets

        size_t Nconst = SignalConstits[counter].size();
        _h["Nconsts"]->fill(Nconst);

        // N-subjettiness -- which one is correct?? the second one does not work...
        // double beta = 1; 
        // for (const PseudoJet& pj : SignalConstits[counter]){
        
        // double tau21, tau32;
        
        // fastjet::contrib::Nsubjettiness nSub1(1, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta));
        // fastjet::contrib::Nsubjettiness nSub2(2, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta));
        // fastjet::contrib::Nsubjettiness nSub3(3, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta));
        
        // double tau1 = nSub1.result(pj);
        // double tau2 = nSub2.result(pj);
        // double tau3 = nSub3.result(pj);
        // (tau1 != 0) ? tau21 = tau2/tau1 : tau21 = -99;
        // (tau2 != 0) ? tau32 = tau3/tau2 : tau32 = -99;
        
        // _h["tau21"]->fill(tau21);
        // _h["tau32"]->fill(tau32);
        // };

        // double beta = 1;        
        // double tau21, tau32;
        
        // fastjet::contrib::Nsubjettiness nSub1(1, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta));
        // fastjet::contrib::Nsubjettiness nSub2(2, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta));
        // fastjet::contrib::Nsubjettiness nSub3(3, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta));
        
        // double tau1 = nSub1.result(j);
        // double tau2 = nSub2.result(j);
        // double tau3 = nSub3.result(j);
        // (tau1 != 0) ? tau21 = tau2/tau1 : tau21 = -99;
        // (tau2 != 0) ? tau32 = tau3/tau2 : tau32 = -99;
        
        // _h["tau21"]->fill(tau21);
        // _h["tau32"]->fill(tau32);
        


        // computes the D values (probabilities)
        _MCbottagger->computeScores(j, tagged_constituents, outputs);

        // distriminant function P for V, H, and top tagger
        double PV=log10(outputs["dnnOutput_V"]/
        (0.9*outputs["dnnOutput_light"]+0.05*outputs["dnnOutput_top"]+0.05*outputs["dnnOutput_H"]));
        double PH=log10(outputs["dnnOutput_H"]/
        (0.9*outputs["dnnOutput_light"]+0.05*outputs["dnnOutput_top"]+0.05*outputs["dnnOutput_V"]));
        double Ptop=log10(outputs["dnnOutput_top"]/
        (0.9*outputs["dnnOutput_light"]+0.05*outputs["dnnOutput_V"]+0.05*outputs["dnnOutput_H"]));

        // fill in the P histograms 
        _h["PV"]->fill(PV);
        _h["PH"]->fill(PH);
        _h["Pt"]->fill(Ptop);

       
       //Tag the jets in an approximation of the MCBot NN.

        MCBot_TagType tag = _MCbottagger->tag(j,tagged_constituents);
        VRCjet_tags.push_back(static_cast<DNN_Category>(tag));

      }


      //Preselection
      //HT > 1250GeV
      double HT = std::accumulate(smeared_small_jets.begin(), smeared_small_jets.end(), 0*GeV, 
                                  [](double count, const Jet& nextjet){return (count+nextjet.pt());});
      if (HT <= 1250*GeV){
        vetoEvent;
      }
      //ETMiss < 200GeV and ETmiss > 40GeV
      if (ETmiss >= 40*GeV || ETmiss <= 200*GeV){
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
      //2 bjets
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
      //Only one of SR and CR is filled (please!)
      string SR = "";
      string VR = "";
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
      
      else if (nVtags + nHtags == 2 && ntoptags >= 2 && smeared_bjets.size()==2){
        _sigBins["XX_2t_2b"]->fill();
      }
      else if (nVtags + nHtags == 2 && ntoptags >= 2 && smeared_bjets.size()>=3){
        _sigBins["XX_2t_3b"]->fill();
      }

      //Now cover the validation regions

      //MultiJet Background Uncertainty Validation Regions (2)
      //TODO: This makes no sense. They say there are two regions but define one?!?
      //Possibly involves altering nVtags? But how?!
      //If I had to guess, by symmetry its HH_1t_2b, but this contradicts the paper. Gong with it for now...
      else if (smeared_bjets.size()==2 && nHtags == 2 && ntoptags == 0 && nVtags == 0){
        _valBins["HH_0t_2b"]->fill();
      }
      else if (smeared_bjets.size()==2 && nHtags == 2 && ntoptags == 1 && nVtags == 0){
        _valBins["HH_1t_2b"]->fill();
      }
      //Closure uncertainty validation regions (7)
      if (nVtags == 2 && nHtags == 0 && ntoptags == 0 && smeared_bjets.size() == 1){
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
      else if (nVtags + nHtags == 2 && ntoptags >= 2 && smeared_bjets.size()==1){
        _valBins["XX_2t_1b"]->fill();
      }

    }
      

      /// Normalise histograms etc., after the run
      void finalize() {
        _h["PV"]->normalize(1);
        _h["PH"]->normalize(1);
        _h["Pt"]->normalize(1);
        _h["LHA"]->normalize(1); 

        // tagger efficiencies and mistag rates...
        
        size_t H_min = _h["PH"]->binIndexAt(0.35);
        size_t V_min = _h["PV"]->binIndexAt(-0.2);
        size_t Top_min = _h["Pt"]->binIndexAt(0.1);

        size_t H_max = _h["PH"]->numBins();
        size_t V_max = _h["PV"]->numBins();
        size_t Top_max = _h["Pt"]->numBins();

        double V = _h["PV"]->integralRange(V_min, V_max-1);
        double H = _h["PH"]->integralRange(H_min, H_max-1);
        double top = _h["Pt"]->integralRange(Top_min, Top_max-1);
        
        // outputs the efficiencies into a file

        std::string OutputFile2 {"Eff_massive_test.csv"};
        std::ofstream file2;
        file2.open(OutputFile2, std::ofstream::app);
        file2 << "V-tag efficiency: " << V << ", " << "H-tag efficiency: " << H << ", " << "Top-tag efficiency: " << top << ", ";
        file2 << "\n";
        file2.close();  

          

      }

      /// @}


      /// @name Histograms and Counters
      /// @{
      map<string, CounterPtr> _sigBins;
      map<string, CounterPtr> _valBins;
      map<string, Histo1DPtr> _h;
      
      
      std::unique_ptr<MCBot_tagger> _MCbottagger;
      size_t _mode;


    };


  RIVET_DECLARE_PLUGIN(ATLAS_2018_I1685207);


}
