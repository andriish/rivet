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


//DEBUG ONLY
template <typename T>
void print_vec(const std::vector<T>& toprint){
  std::cout << "(";
  for (size_t i = 0; i < toprint.size() - 1; ++i){
    std::cout << toprint[i] << ", ";
  }
  std::cout << toprint[toprint.size() - 1] << ")" << std::endl;
}

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



//Very rough and hacky
//Find all Vector bosons, higgs, and tops that are not children of themselves.
//TODO: implement as a projection?
void getVHandtop_fromEvent(std::vector<Particle>& ps, const Event& ge, std::vector<int> wanted_pids = {6,23,24,25}){
  for (const auto& p : ge.allParticles()){
    int pid = p.pid();
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

// angular smearing function: building on JET_SMEAR_ATLAS_RUN2 
Jet JET_SMEAR_ANGULAR(const Jet& j) {
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
  const double fsmear = max(randnorm(1., resolution), 0.); 
  const double mass = j.mass2() > 0 ? j.mass() : 0; //< numerical carefulness...
  
  Jet j1(FourMomentum::mkXYZM(j.px()*fsmear, j.py()*fsmear, j.pz()*fsmear, mass));

  // smearing in eta-phi -- customize the standard deviation in randnorm...
  double dsmear = max(randnorm(0., 0.1), 0.);
  double theta = rand01() * M_2_PI;
  
  return Jet(FourMomentum::mkEtaPhiME(j.eta()+dsmear*cos(theta), mapAngle0To2Pi(j.phi()+dsmear*sin(theta)), j1.mass(), j1.E()));  
}

// angular smearing function: building on JET_SMEAR_ATLAS_RUN2 
// Smears far too much as a sanity check
Jet JET_SMEAR_ANGULAR_ABSURD(const Jet& j) {
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
  const double fsmear = max(randnorm(1., resolution), 0.)*3.5; 
  const double mass = j.mass2() > 0 ? j.mass() : 0; //< numerical carefulness...
  
  Jet j1(FourMomentum::mkXYZM(j.px()*fsmear, j.py()*fsmear, j.pz()*fsmear, mass));

  // smearing in eta-phi -- customize the standard deviation in randnorm...
  double dsmear = max(randnorm(0., 0.1), 0.);
  double theta = rand01() * M_2_PI;
  
  return Jet(FourMomentum::mkEtaPhiME(j.eta()+dsmear*cos(theta), mapAngle0To2Pi(j.phi()+dsmear*sin(theta)), j1.mass(), j1.E()));  
}

/// ATLAS Run 1 jet smearing
  inline Jet JET_SMEAR_ATLAS_RUN1_MAD(const Jet& j) {
    // Jet energy resolution lookup
    //   Implemented by Matthias Danninger for GAMBIT, based roughly on
    //   https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2015-017/
    //   Parameterisation can be still improved, but eta dependence is minimal
    /// @todo Also need a JES uncertainty component?
    static const vector<double> binedges_pt = {0., 50., 70., 100., 150., 200., 1000., 10000.};
    static const vector<double> jer = {0.145, 0.115, 0.095, 0.075, 0.07, 0.05, 0.04, 0.04}; //< note overflow value
    const int ipt = binIndex(j.pT()/GeV, binedges_pt, true);
    if (ipt < 0) return j;
    const double resolution = jer.at(ipt);

    // Smear by a Gaussian centered on 1 with width given by the (fractional) resolution
    /// @todo Is this the best way to smear? Should we preserve the energy, or pT, or direction?
    const double fsmear = max(randnorm(1., resolution), 0.)*3;
    const double mass = j.mass2() > 0 ? j.mass() : 0; //< numerical carefulness...
    Jet rtn(FourMomentum::mkXYZM(j.px()*fsmear, j.py()*fsmear, j.pz()*fsmear, mass));
    //if (deltaPhi(j, rtn) > 0.01) cout << "jdphi: " << deltaPhi(j, rtn) << endl;
    return rtn;
  }


// angular smearing function: that does nothing
Jet JET_SMEAR_FAKE(const Jet& j) {  
  return j;  
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

      // BCG = background mode, SGN = signal mode (V, H, top)  
      // default to signal mode

      _mode = 1; 
      if (getOption("MODE") == "BCG"){
        std::cout << "Mode = BCG" << std::endl;
        _mode = 0;
      }
      else if (getOption("MODE") == "SGN" || getOption("MODE") == "Higgs" ||
              getOption("MODE") == "top" || getOption("MODE") == "Vector") {
        _mode = 1;
        if (getOption("MODE") == "Higgs") {
          _targetParticle = "Higgs";
          std::cout << "Mode = SGN (Higgs)" << std::endl;
        }
        else if (getOption("MODE") == "Vector") {
          _targetParticle = "Vector";
          std::cout << "Mode = SGN (Vector)" << std::endl;
        }
        else if (getOption("MODE") == "top") {
          _targetParticle = "top";
          std::cout << "Mode = SGN (top)" << std::endl;
        }
        else {
          std::cout << "Mode = SGN" << std::endl;
        }
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

      // Histograms for distriminant functions
      book(_h["PV"], "PV", 45,-3,1.5);
      book(_h["PH"], "PH", 55,-3,2.5);
      book(_h["Pt"], "Pt", 55,-3,2.5);

      book(_h["Vt_discriminant"], "Vt_discriminant", 30, -1.5, 1.5);
      book(_h["VH_discriminant"], "VH_discriminant", 30, -1.5, 1.5);
      book(_h["Ht_discriminant"], "Ht_discriminant", 30, -1.5, 1.5);
      book(_h["TripleDiscriminant"], "TripleDiscriminant", 30, -1.5, 1.5);

      book(_h["SigJet_pT"], "SigJet_pT", 15, 40, 2000 );
      book(_h["SigJet_pT_normed"], "SigJet_pT_normed", 15, 40, 2000 );
      book(_s["PV_reweighted"], "PV_reweighted");
      book(_s["PH_reweighted"], "PH_reweighted");
      book(_s["Pt_reweighted"], "Pt_reweighted");
      book(_s["PV_bkgreweighted"], "PV_bkgreweighted");
      book(_s["PH_bkgreweighted"], "PH_bkgreweighted");
      book(_s["Pt_bkgreweighted"], "Pt_bkgreweighted");

      //Jet Mass histograms.
      book(_h["mVRC"], "mVRC", 60, 40, 250);
      book(_h["mVRC_trueHiggs"], "mVRC_trueHiggs", 60, 40, 250);
      book(_h["mVRC_trueVector"], "mVRC_trueVector", 60, 40, 250);
      book(_h["mVRC_trueTop"], "mVRC_trueTop", 60, 40, 250);
      book(_h["mVRC_truebkg"], "mVRC_truebkg", 60, 40, 250);

      book(_h["mVRC_trueHiggsTaggedHiggs"], "mVRC_trueHiggsTaggedHiggs", 60, 40, 250);
      book(_h["mVRC_trueVectorTaggedVector"], "mVRC_trueVectorTaggedVector", 60, 40, 250);
      book(_h["mVRC_trueTopTaggedTop"], "mVRC_trueTopTaggedTop", 60, 40, 250);

      book(_h["mVRC_truebkgTaggedHiggs"], "mVRC_truebkgTaggedHiggs", 60, 40, 250);
      book(_h["mVRC_truebkgTaggedVector"], "mVRC_truebkgTaggedVector", 60, 40, 250);
      book(_h["mVRC_truebkgTaggedTop"], "mVRC_truebkgTaggedTop", 60, 40, 250);
      book(_h["mVRC_truebkgTaggedbkg"], "mVRC_truebkgTaggedbkg", 60, 40, 250);

      //Some interesting things to plot for debugging:
      book(_h["VRC_numConstits"], "VRC_numConstits",{-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5});
      book(_h["sjet1_pt"], "sjet1_pt", 100, -0.0, 2000);
      book(_h["sjet2_pt"], "sjet2_pt", 100, -0.0, 2000);
      book(_h["sjet3_pt"], "sjet3_pt", 100, -0.0, 2000);
      book(_h["sjet1_eta"], "sjet1_eta", 30, -5, 5);
      book(_h["sjet2_eta"], "sjet2_eta", 30, -5, 5);
      book(_h["sjet3_eta"], "sjet3_eta", 30, -5, 5);
      book(_h["sjet1_E"], "sjet1_E", 100, -0.0, 2000);
      book(_h["sjet2_E"], "sjet2_E", 100, -0.0, 2000);
      book(_h["sjet3_E"], "sjet3_E", 100, -0.0, 2000);

      
      //2D histos for post-hoc reweighting
      book(_h2["PVpt"], "PVpt", 45, -3, 1.5, 15, 40, 2000);
      book(_h2["PHpt"], "PHpt", 55, -3, 2.5, 15, 40, 2000);
      book(_h2["Ptpt"], "Ptpt", 55, -3, 2.5, 15, 40, 2000);
      

      // Les Houches angularity
      book(_h["LHA"],"LHA",200, 0, 1);

      //Find the json file
      const std::string nn_datafilename = "ATLAS_2018_I1685207.nn.json.yoda";
      //TODO: Would be nice to use the proper find syntax but there seems to be assumptions
      // about .yoda endings. Someone who understands the paths system better would do it more
      // elegantly.
      const std::string nn_datafilepath =  getDataPath()+"/Rivet/"+nn_datafilename;
      std::cout << "Loading " << nn_datafilepath << std::endl;
      
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
      
      if (_mode == 0){
        signal = FilteredVRCjets;
        SignalConstits = FilteredNewConstits;
      }
      
      else if (_mode  == 1) {

        // get the V, H, and tops from the event
        std::vector<Particle> signalParticles;
        std::vector<Particle> rejectParticles;
        if (_targetParticle == "Higgs"){
          getVHandtop_fromEvent(signalParticles, event, {25}); 
          getVHandtop_fromEvent(rejectParticles, event, {6, 23, 24}); 
        }
        else if (_targetParticle == "Vector"){
          getVHandtop_fromEvent(signalParticles, event, {23, 24}); 
          getVHandtop_fromEvent(rejectParticles, event, {6, 25}); 
        }
        else if (_targetParticle == "top"){
          getVHandtop_fromEvent(signalParticles, event, {6});
          //getVHandtop_fromEvent(rejectParticles, event, {23, 24, 25}); 
        }
        else{
          getVHandtop_fromEvent(signalParticles, event);
        }

        for (size_t i = 0; i < FilteredVRCjets.size(); ++i){
          auto iterator = std::find_if(signalParticles.begin(), signalParticles.end(),
          //TODO: is there a preferred syntax for writing long ugly lambdas?
                          [&FilteredVRCjets, i](const Particle& VHTop){
                            //return (deltaR(momentum3(VHTop.pseudojet()), momentum3(FilteredVRCjets[i])) < 0.1);
                            return (deltaR(momentum3(VHTop.pseudojet()), momentum3(FilteredVRCjets[i])) < 0.75*315/(FilteredVRCjets[i].pt()));
                            });
          if (iterator != signalParticles.end()){ 
            //We don't want jets with multiple tagged particles either of the correct or incorrect type in the jet:
            auto iterator2 = std::find_if(iterator+1, signalParticles.end(),
            //TODO: is there a preferred syntax for writing long ugly lambdas?
                          [&FilteredVRCjets, i](const Particle& VHTop){
                            //return (deltaR(momentum3(VHTop.pseudojet()), momentum3(FilteredVRCjets[i])) < 0.1);
                            return (deltaR(momentum3(VHTop.pseudojet()), momentum3(FilteredVRCjets[i])) < 0.75*315/(FilteredVRCjets[i].pt()));
                            });
            auto iterator3 = std::find_if(rejectParticles.begin(), rejectParticles.end(),
            //TODO: is there a preferred syntax for writing long ugly lambdas?
                          [&FilteredVRCjets, i](const Particle& VHTop){
                            //return (deltaR(momentum3(VHTop.pseudojet()), momentum3(FilteredVRCjets[i])) < 0.1);
                            return (deltaR(momentum3(VHTop.pseudojet()), momentum3(FilteredVRCjets[i])) < 0.75*315/(FilteredVRCjets[i].pt()));
                            });
            
            //Also only consider events where the jet pT in the 40GeV < pt < 2TeV range.
            if (FilteredVRCjets[i].pt() > 40*GeV && FilteredVRCjets[i].pt() < 2000*GeV &&
                 iterator2 == signalParticles.end() && iterator3 == rejectParticles.end()){
              signal.push_back(FilteredVRCjets[i]);
              SignalConstits.push_back(FilteredNewConstits[i]);
            }
          }
        }

        if (signal.size() == 0){
          vetoEvent;
        }

      }

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

        //We need to make sure that the tagged constituents are still pT ordered
        //tagged_constituents.
        isortByPt(tagged_constituents);


        // Les Houches Angularity
        
        double LHA_sum = 0;
        
        for (const PseudoJet& pj : SignalConstits[counter]){
          double pT = pj.pt();
          double dR = pj.squared_distance(j);
          LHA_sum += (pT * pow(dR, 0.25)); 
        }
        
        double LHA = LHA_sum/j.pt();
        
        
        
        //computes the D values (probabilities)
        _MCbottagger->computeScores(j, tagged_constituents, outputs);

        //distriminant function P for V, H, and top tagger
        double PV=log10(outputs["dnnOutput_V"]/
        (0.9*outputs["dnnOutput_light"]+0.05*outputs["dnnOutput_top"]+0.05*outputs["dnnOutput_H"]));
        double PH=log10(outputs["dnnOutput_H"]/
        (0.9*outputs["dnnOutput_light"]+0.05*outputs["dnnOutput_top"]+0.05*outputs["dnnOutput_V"]));
        double Ptop=log10(outputs["dnnOutput_top"]/
        (0.9*outputs["dnnOutput_light"]+0.05*outputs["dnnOutput_V"]+0.05*outputs["dnnOutput_H"]));

        // if (PV < -1.3 || PV > -1.1){
        //   vetoEvent;
        // }

        _h["LHA"]->fill(LHA);

        // fill in the P histograms 
        _h["PV"]->fill(PV);
        _h["PH"]->fill(PH);
        _h["Pt"]->fill(Ptop);

        _h["SigJet_pT"]->fill(j.pt());
        _h["SigJet_pT_normed"]->fill(j.pt());
        

        _h2["PVpt"]->fill(PV, j.pt());
        _h2["PHpt"]->fill(PH, j.pt());
        _h2["Ptpt"]->fill(Ptop, j.pt());

        //Fill in some extra debug info:
        _h["VRC_numConstits"]->fill(tagged_constituents.size());
        _h["sjet1_pt"]->fill(tagged_constituents[0].pt());
        _h["sjet2_pt"]->fill(tagged_constituents.size() > 1 ? tagged_constituents[1].pt() : 0.);
        _h["sjet3_pt"]->fill(tagged_constituents.size() > 2 ? tagged_constituents[2].pt() : 0.);
        _h["sjet1_E"]->fill(tagged_constituents[0].E());
        _h["sjet2_E"]->fill(tagged_constituents.size() > 1 ? tagged_constituents[1].E() : 0.);
        _h["sjet3_E"]->fill(tagged_constituents.size() > 2 ? tagged_constituents[2].E() : 0.);
        _h["sjet1_eta"]->fill(tagged_constituents[0].eta());
        _h["sjet2_eta"]->fill(tagged_constituents.size() > 1 ? tagged_constituents[1].eta() : j.eta());
        _h["sjet3_eta"]->fill(tagged_constituents.size() > 2 ? tagged_constituents[2].eta() : j.eta());

        //Tag the jets in an approximation of the MCBot NN.
        MCBot_TagType DNNtag = _MCbottagger->tag(j,tagged_constituents);
        VRCjet_tags.push_back(static_cast<DNN_Category>(DNNtag));

        //Plot the mass distribution
        _h["mVRC"]->fill(j.m());
        if (_targetParticle == "Higgs"){
          _h["mVRC_trueHiggs"]->fill(j.m());
          if (DNNtag == MCBot_TagType::H){
            _h["mVRC_trueHiggsTaggedHiggs"]->fill(j.m());
          }
        } else if (_targetParticle == "Vector"){
          _h["mVRC_trueVector"]->fill(j.m());
          if (DNNtag == MCBot_TagType::V){
            _h["mVRC_trueVectorTaggedVector"]->fill(j.m());
          }
        }
        else if (_targetParticle == "top"){
          _h["mVRC_trueTop"]->fill(j.m());
          if (DNNtag == MCBot_TagType::top){
            _h["mVRC_trueTopTaggedTop"]->fill(j.m());
          }
        }
        else if (_mode == 0){
          _h["mVRC_truebkg"]->fill(j.m());
          if (DNNtag == MCBot_TagType::H){
            _h["mVRC_truebkgTaggedHiggs"]->fill(j.m());
          }
          else if (DNNtag == MCBot_TagType::V){
            _h["mVRC_truebkgTaggedVector"]->fill(j.m());
          }
          else if (DNNtag == MCBot_TagType::top){
            _h["mVRC_truebkgTaggedTop"]->fill(j.m());
          }
          else {
            _h["mVRC_truebkgTaggedbkg"]->fill(j.m());
          }
        }
        //Plot discriminant function for multi-tag cases, only in the correct signal mode
        //Case 1: is V and H tagged, not T tagged
        if (PV > _MCbottagger->_threshold["PV"] && PH > _MCbottagger->_threshold["PH"] && Ptop < _MCbottagger->_threshold["Pt"]
              && (_targetParticle=="Higgs" || _targetParticle=="Vector")){
          _h["VH_discriminant"]->fill(log10(outputs["dnnOutput_V"]/outputs["dnnOutput_H"]));
        }
        //Case 2: is V and t tagged, not H tagged
        if (PV > _MCbottagger->_threshold["PV"] && PH < _MCbottagger->_threshold["PH"] && Ptop > _MCbottagger->_threshold["Pt"]
              && (_targetParticle=="top" || _targetParticle=="Vector")){
          _h["Vt_discriminant"]->fill(log10(outputs["dnnOutput_V"]/outputs["dnnOutput_top"]));
        }
        //Case 3: is H and t tagged, not H tagged
        if (PV < _MCbottagger->_threshold["PV"] && PH > _MCbottagger->_threshold["PH"] && Ptop > _MCbottagger->_threshold["Pt"]
              && (_targetParticle=="Higgs" || _targetParticle=="top")){
          _h["Ht_discriminant"]->fill(log10(outputs["dnnOutput_H"]/outputs["dnnOutput_top"]));
        }
        //Case 3: triple tagged
        if (PV > _MCbottagger->_threshold["PV"] && PH > _MCbottagger->_threshold["PH"] && Ptop > _MCbottagger->_threshold["Pt"]){
          _h["TripleDiscriminant"]->fill(log10(outputs["dnnOutput_V"]/outputs["dnnOutput_top"]));
        }
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
        if (_h["PV"]->integral() > 0){
          _h["PV"]->normalize(1);
        } if (_h["PH"]->integral() > 0){
          _h["PH"]->normalize(1);
        } if (_h["Pt"]->integral() > 0){
          _h["Pt"]->normalize(1);
        } if (_h["LHA"]->integral() > 0){
          _h["LHA"]->normalize(1);
        }
        if (_h["VH_discriminant"]->integral() > 0){
          _h["VH_discriminant"]->normalize(1);
        } if (_h["Vt_discriminant"]->integral() > 0){
          _h["Vt_discriminant"]->normalize(1);
        } if (_h["Ht_discriminant"]->integral() > 0){
          _h["Ht_discriminant"]->normalize(1);
        } if (_h["TripleDiscriminant"]->integral() > 0){
          _h["TripleDiscriminant"]->normalize(1);
        }
        if (_h["mVRC"]->integral() > 0){
          _h["mVRC"]->normalize(1);
        } if (_h["mVRC_trueHiggs"]->integral() > 0){
          _h["mVRC_trueHiggs"]->normalize(1);
        } if (_h["mVRC_trueVector"]->integral() > 0){
          _h["mVRC_trueVector"]->normalize(1);
        } if (_h["mVRC_trueTop"]->integral() > 0){
          _h["mVRC_trueTop"]->normalize(1);
        } if (_h["mVRC_truebkg"]->integral() > 0){
          _h["mVRC_truebkg"]->normalize(1);
        }
        if (_h["mVRC_trueHiggsTaggedHiggs"]->integral() > 0){
          _h["mVRC_trueHiggsTaggedHiggs"]->normalize(1);
        } if (_h["mVRC_trueVectorTaggedVector"]->integral() > 0){
          _h["mVRC_trueVectorTaggedVector"]->normalize(1);
        } if (_h["mVRC_trueTopTaggedTop"]->integral() > 0){
          _h["mVRC_trueTopTaggedTop"]->normalize(1);
        }
        if (_h["mVRC_truebkgTaggedHiggs"]->integral() > 0){
          _h["mVRC_truebkgTaggedHiggs"]->normalize(1);
        } if (_h["mVRC_truebkgTaggedVector"]->integral() > 0){
          _h["mVRC_truebkgTaggedVector"]->normalize(1);
        }if (_h["mVRC_truebkgTaggedTop"]->integral() > 0){
          _h["mVRC_truebkgTaggedTop"]->normalize(1);
        } if (_h["mVRC_truebkgTaggedbkg"]->integral() > 0){
          _h["mVRC_truebkgTaggedbkg"]->normalize(1);
        }
        if (_h["VRC_numConstits"]->integral() > 0){
          _h["VRC_numConstits"]->normalize(1);
        } if (_h["sjet1_pt"]->integral() > 0){
          _h["sjet1_pt"]->normalize(1);
        } if (_h["sjet2_pt"]->integral() > 0){
          _h["sjet2_pt"]->normalize(1);
        } if (_h["sjet3_pt"]->integral() > 0){
          _h["sjet3_pt"]->normalize(1);
        } if (_h["sjet1_E"]->integral() > 0){
          _h["sjet1_E"]->normalize(1);
        } if (_h["sjet2_E"]->integral() > 0){
          _h["sjet2_E"]->normalize(1);
        } if (_h["sjet3_E"]->integral() > 0){
          _h["sjet3_E"]->normalize(1);
        } if (_h["sjet1_eta"]->integral() > 0){
          _h["sjet1_eta"]->normalize(1);
        } if (_h["sjet2_eta"]->integral() > 0){
          _h["sjet2_eta"]->normalize(1);
        } if (_h["sjet3_eta"]->integral() > 0){
          _h["sjet3_eta"]->normalize(1);
        }
        if (_h["SigJet_pT_normed"]->integral() > 0){
          _h["SigJet_pT_normed"]->normalize(1);
        }





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

        std::string OutputFile2 {"Efficiencies.csv"};
        std::ofstream file2;
        file2.open(OutputFile2, std::ofstream::app);
        file2 << "V-tag efficiency: " << V << ", " << "H-tag efficiency: " << H << ", " << "Top-tag efficiency: " << top << ", ";
        std::cout << "V-tag efficiency: " << V << ", " << "H-tag efficiency: " << H << ", " << "Top-tag efficiency: " << top << ", ";
        file2 << "\n";
        file2.close();  

        if (_h2["PVpt"]->integral() > 0){
          _h2["PVpt"]->normalize(1);
        }
        if (_h2["PHpt"]->integral() > 0){
          _h2["PHpt"]->normalize(1);
        }
        if (_h2["PHpt"]->integral() > 0){
          _h2["Ptpt"]->normalize(1);
        }
        

        //Let's try out the post-hoc reweighting by pT:
        for (const std::string Tag: {"PHpt", "PVpt", "Ptpt"}){
          std::vector<double> PTagEdges =  _h2[Tag]-> xEdges();
          std::vector<double> ptBinEdges = _h2[Tag]-> yEdges();

          
          std::vector<double> PTagCentres = {};
          for (size_t i = 0; i < PTagEdges.size() - 1; ++i){
            PTagCentres.push_back(0.5*(PTagEdges[i] + PTagEdges[i+1]));
          }
          std::vector<double> ptBinCentres = {};
          for (size_t i = 0; i < ptBinEdges.size() - 1; ++i){
            ptBinCentres.push_back(0.5*(ptBinEdges[i] + ptBinEdges[i+1]));
          }


          std::vector<double> pT_volumes = {};
          for (const double pT_value : ptBinCentres){
            double pT_volume = 0;
            for (const double Tag_value : PTagCentres){
              pT_volume += _h2[Tag]->binAt(Tag_value, pT_value).volume();
            }
            pT_volumes.push_back(pT_volume);
          }


          for (size_t i = 0; i < PTagCentres.size(); ++i){
            double Area = 0;
            double BkgReweightArea = 0;
            for (size_t j = 0; j < ptBinCentres.size(); ++j){
              Area += _h2[Tag]->binAt(PTagCentres[i], ptBinCentres[j]).volume() * ( abs(pT_volumes[j]) > 1e-15 ? 1./(pT_volumes[j]) : 0 );
              BkgReweightArea += _h2[Tag]->binAt(PTagCentres[i], ptBinCentres[j]).volume() * ( abs(pT_volumes[j]) > 1e-15 ? _trueSigPtDistrib[j]/(pT_volumes[j]) : 0 );
            }
            //_s[string({Tag[0], Tag[1]})+"_reweighted"]->addPoint(PTagCentres[i], Area, PTagCentres[i] - PTagEdges[i], 0);
            // std::cout << string({Tag[0], Tag[1]})+"_reweighted" << std::endl;
            // std::cout << PTagCentres[i] << ", " << Area << std::endl;
            _s[string({Tag[0], Tag[1]})+"_reweighted"]->addPoint(PTagCentres[i], Area, 5e-2, 0.0001);
            _s[string({Tag[0], Tag[1]})+"_bkgreweighted"]->addPoint(PTagCentres[i], BkgReweightArea, 5e-2, 0.0001);
          }

        }
        

      }

      /// @}


      /// @name Histograms and Counters
      /// @{
      map<string, CounterPtr> _sigBins;
      map<string, CounterPtr> _valBins;
      map<string, Histo1DPtr> _h;
      map<string, Scatter2DPtr> _s;
      map<string, Histo2DPtr> _h2;
      
      
      std::unique_ptr<MCBot_tagger> _MCbottagger;
      size_t _mode;
      std::string _targetParticle;
      std::vector<double> _trueSigPtDistrib {1.826358e-04, 1.524104e-03, 2.572164e-03, 2.158310e-03, 6.498834e-04, 3.057727e-04, 1.395896e-04, 6.725965e-05, 2.747298e-05, 1.371062e-05, 6.156845e-06, 3.673411e-06, 1.241717e-06, 8.795492e-07, 2.069528e-07};


    };


  RIVET_DECLARE_PLUGIN(ATLAS_2018_I1685207);


}
