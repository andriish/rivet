// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/SmearedParticles.hh"
#include "Rivet/Projections/SmearedJets.hh"
#include "Rivet/Tools/RivetLWTNN.hh"

#include "fastjet/contrib/VariableRPlugin.hh"
#include "fastjet/tools/Filter.hh"




namespace Rivet {

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
  Jet JET_SMEAR_ANGULAR_ENERGY_PRESERVED(const Jet& j) {
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
    const double fsmear2 = fsmear*fsmear;
    //Ensure that the smearing doesn't accidently make the mass negative.
    const double newE = ((j.E()*j.E() > j.px()*j.px()*fsmear2 + j.py()*j.py()*fsmear2 + j.pz()*j.pz()*fsmear2) 
                        ? j.E() : sqrt(j.px()*j.px()*fsmear2 + j.py()*j.py()*fsmear2 + j.pz()*j.pz()*fsmear2 + DBL_EPSILON));  
    

    Jet j1(FourMomentum::mkXYZE(j.px()*fsmear, j.py()*fsmear, j.pz()*fsmear, newE));

    // smearing in eta-phi -- customize the standard deviation in randnorm...
    double dsmear = max(randnorm(0., 0.1), 0.);
    double theta = rand01() * M_2_PI;
    
    return Jet(FourMomentum::mkEtaPhiME(j.eta()+dsmear*cos(theta), mapAngle0To2Pi(j.phi()+dsmear*sin(theta)), j1.mass(), j1.E()));  
  }

  // angular smearing function: building on JET_SMEAR_ATLAS_RUN2 
  Jet JET_SMEAR_ANGULAR_PT_PRESERVED(const Jet& j) {
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
    /// @todo TP: This is really guesswork.
    const double fsmear = max(randnorm(1., resolution), 0.); 
    const double fsmear2 = fsmear*fsmear;
    //Ensure that the smearing doesn't accidently make the mass negative.
    //const double newpx = min(j.px()*fsmear, j.pt());
    //const double newpy = sqrt(j.pt2() - newpx*newpx);
    const double newpx = j.px();
    const double newpy = j.py();
    const double newE = (j.E2()*fsmear2 > j.pt2() + j.pz2()*fsmear2) ? j.E() * fsmear : sqrt(j.pt2() + j.pz2() + DBL_EPSILON);
    

    Jet j1(FourMomentum::mkXYZE(newpx, newpy, j.pz()*fsmear, newE));
    // cout << "j: " << j << endl;
    // cout << "j1: " << j1 << endl;

    // smearing in eta-phi -- customize the standard deviation in randnorm...
    double dsmear = max(randnorm(0., 0.1), 0.);
    double theta = rand01() * M_2_PI;
    
    return Jet(FourMomentum::mkEtaPhiME(j.eta()+dsmear*cos(theta), mapAngle0To2Pi(j.phi()+dsmear*sin(theta)), j1.mass(), j1.E()));  
  }

  Jet JET_SMEAR_COMBO(const Jet& j){
    return JET_SMEAR_ATLAS_RUN2(JET_SMEAR_ANGULAR(j));
  }

  enum MCBot_TagType{
    Bkg,
    Vec,
    Higgs,
    top
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


      //Get the analysis mode:
      //STD (default) - 0 -> run analysis normally;
      //VAL-BKG - 1 -> produce background neural network validation plots.
      //VAL-VECTOR - 2 -> produce VECTOR signal neural network validation plots.
      //VAL-HIGGS - 3 -> produce higgs signal neural network validation plots.
      //VAL-TOP - 4 -> produce top signal neural network validation plots.
      //EFF-VAL-BKG - -1 -> Use efficiencies not NN's and produce background neural network validation plots.
      //EFF-VAL-VECTOR - -2 -> Use efficiencies not NN's and produce VECTOR signal neural network validation plots.
      //EFF-VAL-HIGGS - -3 -> Use efficiencies not NN's and produce higgs signal neural network validation plots.
      //EFF-VAL-TOP - -4 -> Use efficiencies not NN's and produce top signal neural network validation plots.
      //EFF-STD - -5 -> Use efficiencies not NN's, but otherwise run analysis normally
      _mode = 0;
      const std::string mode = getOption("MODE");
      if (mode == "VAL-BKG"){
        _mode = 1;
        MSG_DEBUG("Analysis ATLAS_2018_I1685207 running in VAL-BKG mode");
      }
      else if (mode == "VAL-VECTOR"){
        _mode = 2;
        MSG_DEBUG("Analysis ATLAS_2018_I1685207 running in VAL-VECTOR mode");
      }
      else if (mode == "VAL-HIGGS"){
        _mode = 3;
        MSG_DEBUG("Analysis ATLAS_2018_I1685207 running in VAL-HIGGS mode");
      }
      else if (mode == "VAL-TOP"){
        _mode = 4;
        MSG_DEBUG("Analysis ATLAS_2018_I1685207 running in VAL-TOP mode");
      }
      else if (mode == "EFF-STD"){
        _mode = -5;
        MSG_DEBUG("Analysis ATLAS_2018_I1685207 running in EFF-STD mode");
      }
      else if (mode == "EFF-VAL-BKG"){
        _mode = -1;
        MSG_DEBUG("Analysis ATLAS_2018_I1685207 running in EFF-VAL-BKG mode");
      }
      else if (mode == "EFF-VAL-VECTOR"){
        _mode = -2;
        MSG_DEBUG("Analysis ATLAS_2018_I1685207 running in EFF-VAL-VECTOR mode");
      }
      else if (mode == "EFF-VAL-HIGGS"){
        _mode = -3;
        MSG_DEBUG("Analysis ATLAS_2018_I1685207 running in EFF-VAL-HIGGS mode");
      }
      else if (mode == "EFF-VAL-TOP"){
        _mode = -4;
        MSG_DEBUG("Analysis ATLAS_2018_I1685207 running in EFF-VAL-TOP mode");
      }

      //Load the neural net json file.
      //Use hardcode hack for now.
      if (_mode >= 0)
        _nn = mkLWTNN("/home/tomek/PHYSICS_INSTALLS/rivet_seven/rivet/analyses/examples/ATLAS_2018_I1685207.nn.json");

      //Declare projections
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
	    SmearedJets SSj(Sj, JET_SMEAR_ANGULAR_PT_PRESERVED, JET_BTAG_PERFECT);
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

      //Cutflow bins
      for (const string &s : _cutflowNames){
        book(_cutflowBins[s], s);
      }
      
      if (_mode > 0){
        //Validation mode NN histograms
        book(_h["PV"], "PV", 43,-2.8,1.5);
        book(_h["PH"], "PH", 53,-2.8,2.5);
        book(_h["Ptop"], "Ptop", 55,-3,2.5);

        book(_h["Vt_discriminant"], "Vt_discriminant", 30, -1.5, 1.5);
        book(_h["VH_discriminant"], "VH_discriminant", 30, -1.5, 1.5);
        book(_h["Ht_discriminant"], "Ht_discriminant", 30, -1.5, 1.5);
        book(_h["TripleDiscriminant"], "TripleDiscriminant", 30, -1.5, 1.5);

        //Efficiencies And Rejections;
        //Working versions, filled in loops:
        book(_h["tagEfficiency_Vector_in"], "tagEfficiency_Vector_in", {150, 300, 500, 700, 1000, 1500, 2000});
        book(_h["tagEfficiency_Higgs_in"], "tagEfficiency_Higgs_in", {150, 300, 500, 700, 1000, 1500, 2000});
        book(_h["tagEfficiency_Top_in"], "tagEfficiency_Top_in", {150, 300, 500, 700, 1000, 1500, 2000});
        book(_h["bkgRejection_Vector_in"], "bkgRejection_Vector_in", {150, 300, 500, 700, 1000, 1500, 2000});
        book(_h["bkgRejection_Higgs_in"], "bkgRejection_Higgs_in", {150, 300, 500, 700, 1000, 1500, 2000});
        book(_h["bkgRejection_Top_in"], "bkgRejection_Top_in", {150, 300, 500, 700, 1000, 1500, 2000});
        book(_h["tagEfficiency_Vector_out"], "tagEfficiency_Vector_out", {150, 300, 500, 700, 1000, 1500, 2000});
        book(_h["tagEfficiency_Higgs_out"], "tagEfficiency_Higgs_out", {150, 300, 500, 700, 1000, 1500, 2000});
        book(_h["tagEfficiency_Top_out"], "tagEfficiency_Top_out", {150, 300, 500, 700, 1000, 1500, 2000});
        book(_h["bkgRejection_Vector_out"], "bkgRejection_Vector_out", {150, 300, 500, 700, 1000, 1500, 2000});
        book(_h["bkgRejection_Higgs_out"], "bkgRejection_Higgs_out", {150, 300, 500, 700, 1000, 1500, 2000});
        book(_h["bkgRejection_Top_out"], "bkgRejection_Top_out", {150, 300, 500, 700, 1000, 1500, 2000});

        //Final versions, filled during finalise
        book(_s["tagEfficiency_Vector"], "tagEfficiency_Vector");
        book(_s["tagEfficiency_Higgs"], "tagEfficiency_Higgs");
        book(_s["tagEfficiency_Top"], "tagEfficiency_Top");
        book(_s["bkgRejection_Vector"], "bkgRejection_Vector");
        book(_s["bkgRejection_Higgs"], "bkgRejection_Higgs");
        book(_s["bkgRejection_Top"], "bkgRejection_Top");
      }

      // Validation mode Jet Mass histograms.
      if (_mode != 0 && _mode != -5){
        book(_h["mVRC"], "mVRC", 100, 40, 240);
        book(_h["mVRC_trueHiggs"], "mVRC_trueHiggs", 100, 40, 240);
        book(_h["mVRC_trueVector"], "mVRC_trueVector", 100, 40, 240);
        book(_h["mVRC_trueTop"], "mVRC_trueTop", 100, 40, 240);
        book(_h["mVRC_truebkg"], "mVRC_truebkg", 100, 40, 240);

        book(_h["mVRC_trueHiggsTaggedHiggs"], "mVRC_trueHiggsTaggedHiggs", 100, 40, 240);
        book(_h["mVRC_trueVectorTaggedVector"], "mVRC_trueVectorTaggedVector", 100, 40, 240);
        book(_h["mVRC_trueTopTaggedTop"], "mVRC_trueTopTaggedTop", 100, 40, 240);

        book(_h["mVRC_truebkgTaggedHiggs"], "mVRC_truebkgTaggedHiggs", 100, 40, 240);
        book(_h["mVRC_truebkgTaggedVector"], "mVRC_truebkgTaggedVector", 100, 40, 240);
        book(_h["mVRC_truebkgTaggedTop"], "mVRC_truebkgTaggedTop", 100, 40, 240);
        book(_h["mVRC_truebkgTaggedbkg"], "mVRC_truebkgTaggedbkg", 100, 40, 240); 
      }
    }

    
    /// Perform the per-event analysis
    void analyze(const Event& event) {

      //Get particles
      Particles smeared_electrons = apply<ParticleFinder>(event, "SmearedElec").particles();
      Particles smeared_muons = apply<ParticleFinder>(event, "SmearedMuon").particles();
      
      //Get Jets et al
      Jets smeared_small_jets = apply<JetAlg>(event, "smearedSjet").jetsByPt();
      idiscard(smeared_small_jets, Cuts::pt <= 25*GeV && Cuts::abseta >= 2.5); 

      //Overlap removal:
      idiscardIfAnyDeltaRLess(smeared_small_jets, smeared_electrons, 0.2);
      idiscardIfAnyDeltaRLess(smeared_small_jets, smeared_muons, 0.2);
      idiscardIfAny(smeared_small_jets, smeared_muons,
          [](const Jet &j, const Particle &mu){
            return ((std::count_if(j.particles().begin(), j.particles().end(), [](const Particle& tr){
              return tr.pt() > 500*MeV;}) < 3)
              && deltaR(j.p3(), mu.p3()) < 0.4);
          });
      idiscardIfAnyDeltaRLess(smeared_electrons, smeared_small_jets, 0.4);
      idiscardIfAny(smeared_muons, smeared_small_jets, [](const Particle& mu, const Jet& j){
        return deltaR(mu.p3(), j.p3()) < min(0.4, 0.04+10/j.pt());
      });



      // Get info about each-jet's mv2c10 score and do efficiency based b-tagging.
      // TODO: Come up with an elegant (or at least not horrifically ugly) way of preserving this info to the end.
      // FOR NOW: Store as a pair with jet 3-momentum, match back up using deltaR at the end.
      vector<pair<ThreeMomentum,int>> jet_mv2c10_bands;
      jet_mv2c10_bands.reserve(smeared_small_jets.size());
      for (Jet &j : smeared_small_jets){
        jet_mv2c10_bands.push_back({j.p3(), dummy_mv2c10_band(j)});
      }
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
      //Note this is super ugly because we need to keep track of jet substructure by hand.
      PseudoJets TrimmedVRCjets;
      vector<PseudoJets> NewConstits;
      naive_RC_trimmer(VRC_jets, TrimmedVRCjets, 0.05, NewConstits);
      PseudoJet tr_pj;

      PseudoJets UnsortedFilteredVRCjets;
      vector<PseudoJets> UnsortedFilteredNewConstits;
      for (size_t i = 0; i < TrimmedVRCjets.size(); ++i){
        tr_pj = TrimmedVRCjets[i];
        if (tr_pj.pt() > 150*GeV && abs(tr_pj.eta()) < 2.5 && tr_pj.m() > 40*GeV){
          UnsortedFilteredVRCjets.push_back(tr_pj);
          UnsortedFilteredNewConstits.push_back(NewConstits[i]);
        }
      }
      //Ensure Jets are still pT ordered after trimming:
      // Note that the Consituents need to go with it, which leads to some ugliness
      // Feels like there should be an STL algroithm to make this less ugly.
      PseudoJets FilteredVRCjets(UnsortedFilteredVRCjets.size());
      vector<PseudoJets> FilteredNewConstits(UnsortedFilteredVRCjets.size());
      {
        //Build a permutation vector to get the order we want.
        vector<size_t> permutation(UnsortedFilteredVRCjets.size());
        std::iota(permutation.begin(), permutation.end(), 0);
        sort(permutation.begin(), permutation.end(),
         [&UnsortedFilteredVRCjets](size_t i, size_t j){
            return UnsortedFilteredVRCjets[i].pt() > UnsortedFilteredVRCjets[j].pt();
        });
        //Note there's a transform in Rivet namespace also.
        std::transform(permutation.begin(), permutation.end(), 
          FilteredVRCjets.begin(), 
            [&UnsortedFilteredVRCjets](size_t i){return UnsortedFilteredVRCjets[i];
          });
        std::transform(permutation.begin(), permutation.end(), 
          FilteredNewConstits.begin(), [&UnsortedFilteredNewConstits](size_t i){return UnsortedFilteredNewConstits[i];});
      }
      int VRCsize = FilteredVRCjets.size();
      //If there are no vRC jets, why bother carrying on.
      if (VRCsize == 0) {
        vetoEvent;
      }

      // signal - vRC jets for further analysis with their corresponding constituent jets in SignalConstits
      PseudoJets signal; 
      vector<PseudoJets> SignalConstits;
      //Standard or background validation mode: everything is signal.
      if (abs(_mode) < 2 || _mode == -5){
        signal = FilteredVRCjets;
        SignalConstits = FilteredNewConstits;
      }
      else {
        //If running vec/higgs/top mode, need match particle from event record
        //n.b. this is ok as we're not doing "physics" with this, its just validation
        //and matches procedure described in data.
        std::vector<Particle> signalParticles;
        std::vector<Particle> rejectParticles;
        //Higgs
        if (abs(_mode) == 3){
          getVHandtop_fromEvent(signalParticles, event, {25}); 
          getVHandtop_fromEvent(rejectParticles, event, {6, 23, 24}); 
        }
        //Vector
        else if (abs(_mode) == 2){
          getVHandtop_fromEvent(signalParticles, event, {23, 24}); 
          getVHandtop_fromEvent(rejectParticles, event, {6, 25}); 
        }
        //Top
        else if (abs(_mode) == 4){
          getVHandtop_fromEvent(signalParticles, event, {6});
          //getVHandtop_fromEvent(rejectParticles, event, {23, 24, 25}); 
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
            
            //Also only consider events where the jet pT in the 150GeV < pt < 2TeV range.
            // And of mass > 40 GeV, abseta < 2.5
            if (FilteredVRCjets[i].pt() > 150*GeV && FilteredVRCjets[i].pt() < 2000*GeV &&
                 FilteredVRCjets[i].m() > 40*GeV && fabs(FilteredVRCjets[i].eta()) < 2.5 &&
                 iterator2 == signalParticles.end() && iterator3 == rejectParticles.end()){
              signal.push_back(FilteredVRCjets[i]);
              SignalConstits.push_back(FilteredNewConstits[i]);
            }
          }
        }
      }

      if (signal.size() == 0){
        vetoEvent;
      }

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

        //Evaluate the Neural Net and the P-scores.
        //Quick test: Let's only consider leading two jets!
        if (_mode >= 0){
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
        // Alternatively, the NN case
        else {
          //TODO: this is inefficient, but 
          //  - a) The analysis will not be run like this typically
          //  - b) I'll come back and improve it if I really have to
          Particles Vectors, Higgses, Tops;
          getVHandtop_fromEvent(Vectors, event, {23, 24});
          getVHandtop_fromEvent(Higgses, event, {25});
          getVHandtop_fromEvent(Tops, event, {6});

          tag = get_efficiency_tag(j,Vectors,Higgses,Tops);
          JetTags.push_back(tag);
        }
        
        
        //If we're in a validation mode, do validation stuff.
        if (_mode != 0 && _mode != -5){
          if (_mode > 0){
            //distriminant function P for V, H, and top tagger
            const double PV=log10(outputs.at("dnnOutput_V")/
              (0.9*outputs.at("dnnOutput_light")+0.05*outputs.at("dnnOutput_top")+0.05*outputs.at("dnnOutput_H")));
            const double PH=log10(outputs.at("dnnOutput_H")/
              (0.9*outputs.at("dnnOutput_light")+0.05*outputs.at("dnnOutput_top")+0.05*outputs.at("dnnOutput_V")));
            const double Ptop=log10(outputs.at("dnnOutput_top")/
              (0.9*outputs.at("dnnOutput_light")+0.05*outputs.at("dnnOutput_V")+0.05*outputs.at("dnnOutput_H")));

            _h["PV"]->fill(PV);
            _h["PH"]->fill(PH);
            _h["Ptop"]->fill(Ptop);

            //If we're in backgraound validation mode, fill background rejection plots
            if (_mode == 1){
              if (tag == MCBot_TagType::Vec)
                _h["bkgRejection_Vector_in"]->fill(j.pt());
              else  
                _h["bkgRejection_Vector_out"]->fill(j.pt());
              
              if (tag == MCBot_TagType::Higgs)
                _h["bkgRejection_Higgs_in"]->fill(j.pt());
              else  
                _h["bkgRejection_Higgs_out"]->fill(j.pt());

              if (tag == MCBot_TagType::top)
                _h["bkgRejection_Top_in"]->fill(j.pt());
              else  
                _h["bkgRejection_Top_out"]->fill(j.pt());
            }


            //Tag Efficiencies & MultiTag discriminant functions (see fig 2 & 4 of paper)
            // Does not apply to bkg jets
            if ( _mode > 1 ){
              //Tag efficiencies
              if (_mode == 2 && tag == MCBot_TagType::Vec)
                _h["tagEfficiency_Vector_in"]->fill(j.pt());
              else if (_mode == 2)
                _h["tagEfficiency_Vector_out"]->fill(j.pt());
              if (_mode == 3 && tag == MCBot_TagType::Higgs)
                _h["tagEfficiency_Higgs_in"]->fill(j.pt());
              else if (_mode == 3)
                _h["tagEfficiency_Higgs_out"]->fill(j.pt());
              if (_mode == 4 && tag == MCBot_TagType::top)
                _h["tagEfficiency_Top_in"]->fill(j.pt());
              else if (_mode == 4)
                _h["tagEfficiency_Top_out"]->fill(j.pt());
              


              // discriminant functions
              //Case 1: The PV function
              if (PV > _threshold.at("PV") && PH > _threshold.at("PH") && Ptop < _threshold.at("Pt")
                && (_mode == 2 || _mode == 3)){
                  _h["VH_discriminant"]->fill(log10(outputs.at("dnnOutput_V")/outputs.at("dnnOutput_H")));
              }
              //Case 2: The Vt function.
              if (PV > _threshold.at("PV") && PH < _threshold.at("PH") && Ptop > _threshold.at("Pt")
                && (_mode == 2 || _mode == 4 )){
                _h["Vt_discriminant"]->fill(log10(outputs.at("dnnOutput_V")/outputs.at("dnnOutput_top")));
              }
              //Case 3: The Ht function.
              if (PV < _threshold.at("PV") && PH > _threshold.at("PH") && Ptop > _threshold.at("Pt")
                && ( _mode == 3 || _mode == 4 )){
                  _h["Ht_discriminant"]->fill(log10(outputs.at("dnnOutput_H")/outputs.at("dnnOutput_top")));
              }
              //Case 4: The triple-tag plot.
              if (PV > _threshold.at("PV") && PH > _threshold.at("PH") && Ptop > _threshold.at("Pt")){
                _h["TripleDiscriminant"]->fill(log10(outputs.at("dnnOutput_V")/outputs.at("dnnOutput_top")));
              }
            }
          }

          //Variable-R Reclustered Jet Masses:
          _h["mVRC"]->fill(j.m());
          if (abs(_mode) == 3){
            _h["mVRC_trueHiggs"]->fill(j.m());
            if (tag == MCBot_TagType::Higgs){
              _h["mVRC_trueHiggsTaggedHiggs"]->fill(j.m());
            }
          } else if(abs(_mode) == 2){
            _h["mVRC_trueVector"]->fill(j.m());
            if (tag == MCBot_TagType::Vec){
              _h["mVRC_trueVectorTaggedVector"]->fill(j.m());
            }
          }
          else if (abs(_mode) == 4){
            _h["mVRC_trueTop"]->fill(j.m());
            if (tag == MCBot_TagType::top){
              _h["mVRC_trueTopTaggedTop"]->fill(j.m());
            }
          }
          else if (abs(_mode) == 1){
            _h["mVRC_truebkg"]->fill(j.m());
            if (tag == MCBot_TagType::Higgs){
              _h["mVRC_truebkgTaggedHiggs"]->fill(j.m());
            }
            else if (tag == MCBot_TagType::Vec){
              _h["mVRC_truebkgTaggedVector"]->fill(j.m());
            }
            else if (tag == MCBot_TagType::top){
              _h["mVRC_truebkgTaggedTop"]->fill(j.m());
            }
            else {
              _h["mVRC_truebkgTaggedbkg"]->fill(j.m());
            }
          }

        }

      }

      

      //Now let's get on with the actual analysis
      //Preselection
      if (smeared_electrons.size() != 0){
        vetoEvent;
      }
      if (smeared_muons.size() != 0){
        vetoEvent;
      }
      //HT > 1250GeV
      // HT = total scalar sum of the transverse momenta of all track particles and energy deposits
      double HT = std::accumulate(smeared_small_jets.begin(), smeared_small_jets.end(), 0*GeV, 
                                  [](double count, const Jet& nextjet){return (count+nextjet.pt());}); 
      if (HT <= 1250*GeV){
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

      //"preselection" cutflow should be.
      _cutflowBins["presel"]->fill();

      // TODO: I believe this cutf is superfluous for the main analysis
      // (i.e. it cuts nothing that isn't cut elsewhere), but for validation
      // we need the cutflow.
      // It's probably a slightly inefficiecnt ordering to leave it in though.
      // Cut for >=2 boson candiates:
      size_t nCandidates = 0;
      for (size_t i = 0; i < smeared_small_jets.size(); ++i){
        if (JetTags[i] == MCBot_TagType::Vec || JetTags[i] == MCBot_TagType::Higgs){
          ++nCandidates;
        } 
        else if (69*GeV < smeared_small_jets[i].mass() && smeared_small_jets[i].mass() < 155*GeV){
          ++nCandidates;
        }
      }
      if (nCandidates < 2){
        vetoEvent;
      }
      _cutflowBins[">=2 boson candidates"]->fill();

      //2 bjets
      if (smeared_bjets.size() < 2){
        vetoEvent;
      }
      _cutflowBins[">=2 b tags"]->fill();

      //ETMiss < 200GeV and ETmiss > 40GeV
      if (ETmiss <= 40*GeV || ETmiss >= 200*GeV){
        vetoEvent;
      }
      //n.b the cutflow table only mentions ETMiss > 40, but assume it means both?
      _cutflowBins["ETmiss >= 40"]->fill();

      // //Two vRC jets tagged V or H
      // int nVtags = std::count_if(JetTags.begin(), JetTags.end(), 
      //                             [](const MCBot_TagType tt){return (tt == MCBot_TagType::Vec);});
      // int nHtags = std::count_if(JetTags.begin(), JetTags.end(), 
      //                             [](const MCBot_TagType tt){return (tt == MCBot_TagType::Higgs);});
      // int ntoptags = std::count_if(JetTags.begin(), JetTags.end(), 
      //                             [](const MCBot_TagType tt){return (tt == MCBot_TagType::top);});

      // Get two highest pT tags (JetTags should be pT ordered as Signal is pT ordered)
      int nVtags = 0, nHtags = 0, ntoptags = 0;
      for (const MCBot_TagType tt : JetTags){
        switch(tt){
          case (MCBot_TagType::Vec):
            ++nVtags; break;
          case (MCBot_TagType::Higgs):
            ++nHtags; break;
          case (MCBot_TagType::top):
            ++ntoptags; break;
          default:
            break; //Do nothing
        }
        if (nVtags + nHtags == 2) break;
      }

      if (nVtags + nHtags < 2){
        vetoEvent;
      }

      _cutflowBins[">=2 bosons (DNN)"]->fill();

      //select signal/validation/control region.
      //TODO: The pre-selection cuts say 2 or more (v or H) tagged jets, but each signal region requires only two.
      //  TODO: I understand the above now, needs a slight rewrite. Should lead to an increase in counts?
      // I also can't see some sort of tie-break procedure outlined (e.g. take tags of two highest pT jets)

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
      //TODO: Reading the int note, it's definitely >= 2. Come back to this later.
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
      //Normalise sig/val region counts to XS
      //TODO: luminosityfb(), luminosity() don't seem to work.
      //Hardcoded for now.
      const double sf = 1000*crossSection()*(36.1)/sumW();
      for (const string& s : _sigRegionNames ){
        _sigBins[s]->scaleW(sf);
      }
      for (const string& s : _controlRegionNames){
        _valBins[s]->scaleW(sf);
      }
      for (const string& s : _cutflowNames){
        _cutflowBins[s]->scaleW(sf);
      }

      //Calculate efficiencies/rejections (NN, Vaidation only):
      // rejections
      if (_mode == 1){
        for (size_t i = 0; i < _h["bkgRejection_Vector_out"]->bins().size(); ++i){
          _s["bkgRejection_Vector"]->addPoint(_h["bkgRejection_Vector_out"]->bins()[i].xMid(), 
                                                  ((_h["bkgRejection_Vector_out"]->bins()[i].sumW()+_h["bkgRejection_Vector_in"]->bins()[i].sumW())
                                                    /_h["bkgRejection_Vector_in"]->bins()[i].sumW()));

          _s["bkgRejection_Higgs"]->addPoint(_h["bkgRejection_Higgs_out"]->bins()[i].xMid(), 
                                                  ((_h["bkgRejection_Higgs_out"]->bins()[i].sumW()+_h["bkgRejection_Higgs_in"]->bins()[i].sumW())
                                                    /_h["bkgRejection_Higgs_in"]->bins()[i].sumW()));

          _s["bkgRejection_Top"]->addPoint(_h["bkgRejection_Top_out"]->bins()[i].xMid(), 
                                                  ((_h["bkgRejection_Top_out"]->bins()[i].sumW()+_h["bkgRejection_Top_in"]->bins()[i].sumW())
                                                    /_h["bkgRejection_Top_in"]->bins()[i].sumW()));
        }
      }
      else if (_mode == 2){
        //Efficiency for vectors
        for (size_t i = 0; i < _h["tagEfficiency_Vector_in"]->bins().size(); ++i){
          _s["tagEfficiency_Vector"]->addPoint(_h["tagEfficiency_Vector_in"]->bins()[i].xMid(), 
                                                  _h["tagEfficiency_Vector_in"]->bins()[i].sumW()/
                                                    (_h["tagEfficiency_Vector_in"]->bins()[i].sumW()+_h["tagEfficiency_Vector_out"]->bins()[i].sumW()));
        }
      }
      else if (_mode == 3){
        //Efficiency for Higgs
        for (size_t i = 0; i < _h["tagEfficiency_Higgs_in"]->bins().size(); ++i){
          _s["tagEfficiency_Higgs"]->addPoint(_h["tagEfficiency_Higgs_in"]->bins()[i].xMid(), 
                                                  _h["tagEfficiency_Higgs_in"]->bins()[i].sumW()/
                                                    (_h["tagEfficiency_Higgs_in"]->bins()[i].sumW()+_h["tagEfficiency_Higgs_out"]->bins()[i].sumW()));
        }
      }
      else if (_mode == 4){
        //Efficiency for Higgs
        for (size_t i = 0; i < _h["tagEfficiency_Top_in"]->bins().size(); ++i){
          _s["tagEfficiency_Top"]->addPoint(_h["tagEfficiency_Top_in"]->bins()[i].xMid(), 
                                                  _h["tagEfficiency_Top_in"]->bins()[i].sumW()/
                                                    (_h["tagEfficiency_Top_in"]->bins()[i].sumW()+_h["tagEfficiency_Top_out"]->bins()[i].sumW()));
        }
      }
      

      //Validation Plot Normalisation.
      if (_mode > 0){
        if (_h["PV"]->integral() > 0){
          _h["PV"]->normalize(1);
        } if (_h["PH"]->integral() > 0){
          _h["PH"]->normalize(1);
        } if (_h["Ptop"]->integral() > 0){
          _h["Ptop"]->normalize(1);
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
      }
      if (_mode != 0 && _mode !=-5) {
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
      }
    }


    /// @}

  
    /// @name methods

    private:
    //trims reclustered jets by removing all subjets with pt < minpt. Does not work in place,
    //has seperate inout&output. Also returns a vector of new constituents so we definitely save this
    // info before the cluster sequence pointer disappears in a puff of smoke.
    //If I can use a fastjet trimmer to do this more elegantly please let me know.
    static void naive_RC_trimmer(const Rivet::PseudoJets& input, Rivet::PseudoJets& output,
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

    MCBot_TagType getTag(const map<string, double> &outputs){

      const double PV=log10(outputs.at("dnnOutput_V")/
        (0.9*outputs.at("dnnOutput_light")+0.05*outputs.at("dnnOutput_top")+0.05*outputs.at("dnnOutput_H")));
      const double PH=log10(outputs.at("dnnOutput_H")/
        (0.9*outputs.at("dnnOutput_light")+0.05*outputs.at("dnnOutput_top")+0.05*outputs.at("dnnOutput_V")));
      const double Ptop=log10(outputs.at("dnnOutput_top")/
        (0.9*outputs.at("dnnOutput_light")+0.05*outputs.at("dnnOutput_V")+0.05*outputs.at("dnnOutput_H")));

        
      const bool hasVtag = PV > _threshold.at("PV");
      const bool hasHtag = PH > _threshold.at("PH");
      const bool hasttag = Ptop > _threshold.at("Pt");

      //No tags => Backgroundjet
      if (!(hasVtag || hasHtag || hasttag))
        return MCBot_TagType::Bkg;

      //Each of the single tag cases (and triple tagged => Higgs)
      if (hasVtag && !(hasHtag || hasttag))
        return MCBot_TagType::Vec;
      if ((hasHtag && !(hasVtag || hasttag)) || (hasHtag && hasVtag && hasttag))
        return MCBot_TagType::Higgs;
      if (hasttag && !(hasVtag || hasHtag))
        return MCBot_TagType::top;

      //Cases with two tags - this requires more work.
      if (hasVtag && hasHtag){
        const double disc_score = log10(outputs.at("dnnOutput_V")/outputs.at("dnnOutput_H"));
        return (disc_score > _tiebreak_thresholds.at("H_V")) ? MCBot_TagType::Vec : MCBot_TagType::Higgs;
      }
      if (hasVtag && hasttag){
        const double disc_score = log10(outputs.at("dnnOutput_V")/outputs.at("dnnOutput_top"));
        return (disc_score > _tiebreak_thresholds.at("t_V")) ? MCBot_TagType::Vec : MCBot_TagType::top;
      }
      if (hasHtag && hasttag){
        const double disc_score = log10(outputs.at("dnnOutput_H")/outputs.at("dnnOutput_top"));
        return (disc_score > _tiebreak_thresholds.at("t_H")) ? MCBot_TagType::Higgs : MCBot_TagType::top;
      }

      //ELSE - SHOULDN'T EVER REACH THIS:
      MSG_ERROR("MCBot tagging not succesful -- debugging required");
      throw Error("MCBot tagging not succesful -- debugging required");
    }

    /// @}

    /// @name Efficiency based taggers (using fig 2)
    static MCBot_TagType get_efficiency_bkg_tag(const double jetpt){
      //TODO: How to interpret the background "rejection" to me is a
      // little ambigious, particularly in the context of double or triple tagging
      double cut_off_pby_vec, cut_off_pby_higgs, cut_off_pby_top;
      if (jetpt < 150){
        cut_off_pby_vec = 0;
        cut_off_pby_higgs = 0;
        cut_off_pby_top = 0;
      } else if (jetpt < 300){
        cut_off_pby_vec = 1/6.3;
        cut_off_pby_higgs = 1/33.;
        cut_off_pby_top = 1/16.5;
      } else if (jetpt < 500){
        cut_off_pby_vec = 1/6.9;
        cut_off_pby_higgs = 1/39.;
        cut_off_pby_top = 1/17.;
      } else if (jetpt < 700){
        cut_off_pby_vec = 1/6.2;
        cut_off_pby_higgs = 1/45.;
        cut_off_pby_top = 1/18.1;
      } else if (jetpt < 1000){
        cut_off_pby_vec = 1/4.8;
        cut_off_pby_higgs = 1/46.;
        cut_off_pby_top = 1/18.75;
      } else if (jetpt < 1500){
        cut_off_pby_vec = 1/4.5;
        cut_off_pby_higgs = 1/43.;
        cut_off_pby_top = 1/18.2;
      } else {
        cut_off_pby_vec = 1/4.9;
        cut_off_pby_higgs = 1/38.;
        cut_off_pby_top = 1/12.8;
      }
      const double vecrndm = rand01();
      const double higgsrndm = rand01();
      const double toprndm = rand01();

      const bool vectag = cut_off_pby_vec < vecrndm;
      const bool higgstag = cut_off_pby_higgs < higgsrndm;
      const bool toptag = cut_off_pby_top < toprndm;

      if (!vectag && !higgstag && !toptag){
         return MCBot_TagType::Bkg;
      }
      else if (vectag && !higgstag && !toptag){
        return MCBot_TagType::Vec;
      }
      else if (!vectag && higgstag && !toptag){
        return MCBot_TagType::Higgs;
      }
      else if (!vectag && !higgstag && toptag){
        return MCBot_TagType::top;
      }
      else if (vectag && higgstag && !toptag){
        return (rand01() > 0.5 ? MCBot_TagType::Vec : MCBot_TagType::Higgs);
      }
      else if (vectag && !higgstag && toptag){
        return (rand01() > 0.5 ? MCBot_TagType::Vec : MCBot_TagType::top);
      }
      else if (!vectag && higgstag && toptag){
        return (rand01() > 0.5 ? MCBot_TagType::Higgs : MCBot_TagType::top);
      }
      else {
        const double score = rand01();
        if (score <  0.333)
          return MCBot_TagType::Vec;
        else if (score < 0.667)
          return MCBot_TagType::Higgs;
        else
          return MCBot_TagType::top;
      }
    }
    //Tag a jet that is "truly" a vector with either a vec or bkg tag
    static MCBot_TagType get_efficiency_vec_tag(const double jetpt){
      double cut_off_pby;
      if (jetpt < 150){
        cut_off_pby = 0;
      } else if (jetpt < 300){
        cut_off_pby = 0.58;
      } else if (jetpt < 500){
        cut_off_pby = 0.7;
      } else if (jetpt < 700){
        cut_off_pby = 0.75;
      } else if (jetpt < 1000){
        cut_off_pby = 0.74;
      } else if (jetpt < 1500){
        cut_off_pby = 0.7;
      } else {
        cut_off_pby = 0.62;
      }

      if (rand01() < cut_off_pby)
        return MCBot_TagType::Vec;
      else return MCBot_TagType::Bkg;
    }
    //Tag a jet that is "truly" a higgs with either a higgs or bkg tag
    static MCBot_TagType get_efficiency_higgs_tag(const double jetpt){
      double cut_off_pby;
      if (jetpt < 150){
        cut_off_pby = 0;
      } else if (jetpt < 300){
        cut_off_pby = 0.68;
      } else if (jetpt < 500){
        cut_off_pby = 0.73;
      } else if (jetpt < 700){
        cut_off_pby = 0.69;
      } else if (jetpt < 1000){
        cut_off_pby = 0.62;
      } else if (jetpt < 1500){
        cut_off_pby = 0.58;
      } else {
        cut_off_pby = 0.65;
      }

      if (rand01() < cut_off_pby)
        return MCBot_TagType::Higgs;
      else return MCBot_TagType::Bkg;
    }
    //Tag a jet that is "truly" a top with either a top or bkg tag
    static MCBot_TagType get_efficiency_top_tag(const double jetpt){
      double cut_off_pby;
      if (jetpt < 150){
        cut_off_pby = 0;
      } else if (jetpt < 300){
        cut_off_pby = 0.53;
      } else if (jetpt < 500){
        cut_off_pby = 0.62;
      } else if (jetpt < 700){
        cut_off_pby = 0.63;
      } else if (jetpt < 1000){
        cut_off_pby = 0.59;
      } else if (jetpt < 1500){
        cut_off_pby = 0.57;
      } else {
        cut_off_pby = 0.56;
      }

      if (rand01() < cut_off_pby)
        return MCBot_TagType::top;
      else return MCBot_TagType::Bkg;
    }


    MCBot_TagType get_efficiency_tag(const Jet &j, const Particles& Vectors, const Particles& Higgses, const Particles& Tops){
      auto vec_iterator = std::find_if(Vectors.begin(), Vectors.end(),
          //TODO: is there a preferred syntax for writing long ugly lambdas?
                          [&j](const Particle& VHTop){
                            return (deltaR(momentum3(VHTop.pseudojet()), momentum3(j)) < 0.75*315*GeV/(j.pt()));
                          });
      auto higgs_iterator = std::find_if(Higgses.begin(), Higgses.end(),
          //TODO: is there a preferred syntax for writing long ugly lambdas?
                          [&j](const Particle& VHTop){
                            return (deltaR(momentum3(VHTop.pseudojet()), momentum3(j)) < 0.75*315*GeV/(j.pt()));
                          });
      auto top_iterator = std::find_if(Tops.begin(), Tops.end(),
          //TODO: is there a preferred syntax for writing long ugly lambdas?
                          [&j](const Particle& VHTop){
                            return (deltaR(momentum3(VHTop.pseudojet()), momentum3(j)) < 0.75*315*GeV/(j.pt()));
                          });

      const bool isVecJet = (vec_iterator != Vectors.end());
      const bool isHiggsJet = (higgs_iterator != Higgses.end());
      const bool isTopJet = (top_iterator != Tops.end());

      if (!isVecJet && !isHiggsJet && !isTopJet)
        return get_efficiency_bkg_tag(j.pt());
      else if (isVecJet && !isHiggsJet && !isTopJet)
        return get_efficiency_vec_tag(j.pt());
      else if (!isVecJet && isHiggsJet && !isTopJet)
        return get_efficiency_higgs_tag(j.pt());
      else if (!isVecJet && !isHiggsJet && isTopJet)
        return get_efficiency_top_tag(j.pt());
      
      //The efficiencies almost certainly break down in these cases
      // But I gotta do something
      // note to self: this is a great demonstration of why giving us 
      // the network is better - remember to put in the presentation.
      else if (isVecJet && isHiggsJet && !isTopJet)
        return (rand01() > 0.5 ? get_efficiency_vec_tag(j.pt()) : get_efficiency_higgs_tag(j.pt()));
      else if (isVecJet && !isHiggsJet && isTopJet)
        return (rand01() > 0.5 ? get_efficiency_vec_tag(j.pt()) : get_efficiency_top_tag(j.pt()));
      else if (!isVecJet && isHiggsJet && isTopJet)
        return (rand01() > 0.5 ? get_efficiency_higgs_tag(j.pt()) : get_efficiency_top_tag(j.pt()));

      //Triple case is a weird jet indeed.
      else return (rand01() > 0.333 ? get_efficiency_vec_tag(j.pt()) : rand01() > 0.5 ? get_efficiency_higgs_tag(j.pt()) : get_efficiency_top_tag(j.pt()));
      
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
      const bool btag  = rand01() < beff77;
      // Remove b-tags if needed, and add a dummy one if needed
      if (!btag && j.bTagged()) j.tags().erase(std::remove_if(j.tags().begin(), j.tags().end(), hasBottom), j.tags().end());
      if (btag && !j.bTagged()) j.tags().push_back(Particle(PID::BQUARK, j.mom())); ///< @todo Or could use the/an actual clustered b-quark momentum?
      
      //Ignore c-tagging.

      if (variate < beff60) return 60;
      else if (variate < beff70) return 70;
      else if (variate < beff77) return 77;
      else if (variate < beff85) return 85;
      else return 100;
    }

    double mv2c10_score_from_band(const int band){
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


    /// @}

    /// @name Member variables
    map<string, CounterPtr> _sigBins;
    map<string, CounterPtr> _valBins;
    map<string, CounterPtr> _cutflowBins;
    map<string, Histo1DPtr> _h;
    map<string, Scatter2DPtr> _s;
    map<string, Histo2DPtr> _h2;

    int _mode;
    std::unique_ptr<lwt::LightweightNeuralNetwork> _nn;  
    const std::map<string, double> _threshold = {{"PV", -0.2}, {"PH", 0.35}, {"Pt", 0.1}};
    const std::map<string, double> _tiebreak_thresholds = {{"t_V", -0.3}, {"H_V", -0.55}, {"t_H", 0.2}};
    //List of sig/control region names for easy iteration:
    const std::vector<string> _sigRegionNames = {"VV_0t_2b", "VV_0t_3b", "VV_1t_2b", "VV_1t_3b",
                                                 "VH_0t_2b", "VH_0t_3b", "VH_1t_2b", "VH_1t_3b",
                                                 "HH_0t_3b", "HH_1t_3b", "XX_2t_2b", "XX_2t_3b"};
    const std::vector<string> _controlRegionNames = {"HH_0t_2b", "HH_1t_2b",
                                                     "VV_0t_1b", "VV_1t_1b",
                                                     "VH_0t_1b", "VH_1t_1b",
                                                     "HH_0t_1b", "HH_1t_1b", "XX_2t_1b"};
    const vector<string> _cutflowNames = {"presel", ">=2 boson candidates", ">=2 b tags", "ETmiss >= 40", ">=2 bosons (DNN)" };

    /// @}

  };

 

  RIVET_DECLARE_PLUGIN(ATLAS_2018_I1685207);

}
