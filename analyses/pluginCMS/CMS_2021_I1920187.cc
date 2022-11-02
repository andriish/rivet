// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Jet.hh"

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Recluster.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include <algorithm>


using std::vector;
using namespace fastjet;

/// @brief Study of quark and gluon jet substructure in Z+jet and dijet events from pp collisions

namespace Rivet {

  /// @brief Routine for QG substructure analysis
  class CMS_2021_I1920187 : public Analysis {
  public:

    /// Constructor
    CMS_2021_I1920187()
      : Analysis("CMS_2021_I1920187")
    {}

    /// Book histograms and initialise projections before the run
    void init() {

      _mode = 0;
      if ( getOption("MODE") == "DIJET" ) _mode = 1;
      else if ( getOption("MODE") == "ZJET" ) _mode = 2;
      else {
	MSG_WARNING("Mode not specified in CMS_2021_I1920187, using DIJET");
	_mode = 1;
      }

      // Initialise and register projections
      FinalState fs(Cuts::abseta < 5 && Cuts::pT > 0*GeV);
      // Z-jet
      if(_mode==2) {
	// for the muons
	double mu_pt = 26.;
	double mz_min = (90-20);
	double mz_max = (90+20);
	double eta_max = 2.4;
	ZFinder zfinder(fs,
			Cuts::pT > mu_pt*GeV  && Cuts::abseta < eta_max,
			PID::MUON,
			mz_min*GeV, mz_max*GeV,
			0.1, ZFinder::ClusterPhotons::NONE, ZFinder::AddPhotons::NO);
	declare(zfinder, "ZFinder");

	eta_max = 2.4;
	FinalState fs_muons(Cuts::abseta < eta_max && Cuts::pT > 0*GeV);
	IdentifiedFinalState muons_noCut(fs_muons, {PID::MUON, PID::ANTIMUON});
	declare(muons_noCut, "MUONS_NOCUT");
	// Particles for the jets
	VetoedFinalState jet_input(fs);
	jet_input.vetoNeutrinos();
	jet_input.addVetoOnThisFinalState(getProjection<ZFinder>("ZFinder"));
	declare(jet_input, "JET_INPUT");
	_ptBinsGen = { 50, 65, 88, 120, 150, 186, 254, 326, 408, 1500};
      }
      // dijet
      else {
	// Particles for the jets
	VetoedFinalState jet_input(fs);
	jet_input.vetoNeutrinos();
	declare(jet_input, "JET_INPUT");
	_ptBinsGen = {50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 1000, 4000};
      }
      // Book histograms
      // resize vectors appropriately
      uint nHistsRadii = _jetRadii.size();
      uint nHistsLambda = _lambdaVars.size();
      uint nHistsPt = _ptBinsGen.size()-1;

      // Z jet
      if (_mode==2) {
	_h_zpj.resize(nHistsRadii, vector<vector<Histo1DPtr> >(nHistsLambda, vector<Histo1DPtr>(nHistsPt)));
	_h_zpj_groomed.resize(nHistsRadii, vector<vector<Histo1DPtr> >(nHistsLambda, vector<Histo1DPtr>(nHistsPt)));

	// Now book histos
	// remember 1-indexed

	// yoda plot naming scheme
	// --------------------------------------------------------------------------
	// channel (ak4/8 radii [000, 100] + zjet / dijet cen / fwd [00, 10, 20] + groomed versions [0, 1])
	// lambda variable; neutral+charged & charged-only are treated separately [0000,..,4000]
	// pT bin [00000,..,120000]
	for (uint radiusInd=0; radiusInd < _jetRadii.size(); radiusInd++) {
	  for (uint lambdaInd=0; lambdaInd < _lambdaVars.size(); lambdaInd++) {
	    for (uint ptInd=0; ptInd < _ptBinsGen.size()-1; ptInd++) {
	      book(_h_zpj[radiusInd][lambdaInd][ptInd], _hepdata_index[0+100*radiusInd+1000*lambdaInd+10000*ptInd], 1, 1);
	      book(_h_zpj_groomed[radiusInd][lambdaInd][ptInd], _hepdata_index[10+100*radiusInd+1000*lambdaInd+10000*ptInd], 1, 1);
	    }
	  }
	}
      }
      // di jet
      else {
	_h_dijet_cen.resize(nHistsRadii, vector<vector<Histo1DPtr> >(nHistsLambda, vector<Histo1DPtr>(nHistsPt)));
	_h_dijet_cen_groomed.resize(nHistsRadii, vector<vector<Histo1DPtr> >(nHistsLambda, vector<Histo1DPtr>(nHistsPt)));
	_h_dijet_fwd.resize(nHistsRadii, vector<vector<Histo1DPtr> >(nHistsLambda, vector<Histo1DPtr>(nHistsPt)));
	_h_dijet_fwd_groomed.resize(nHistsRadii, vector<vector<Histo1DPtr> >(nHistsLambda, vector<Histo1DPtr>(nHistsPt)));

	// Now book histos
	// remember 1-indexed

	// yoda plot naming scheme
	// --------------------------------------------------------------------------
	// channel (ak4/8 radii [000, 100] + zjet / dijet cen / fwd [00, 10, 20] + groomed versions [0, 1])
	// lambda variable; neutral+charged & charged-only are treated separately [0000,..,4000]
	// pT bin [00000,..,120000]
	for (uint radiusInd=0; radiusInd < _jetRadii.size(); radiusInd++) {
	  for (uint lambdaInd=0; lambdaInd < _lambdaVars.size(); lambdaInd++) {
	    for (uint ptInd=0; ptInd < _ptBinsGen.size()-1; ptInd++) {
	      book(_h_dijet_cen[radiusInd][lambdaInd][ptInd], _hepdata_index[1+100*radiusInd+1000*lambdaInd+10000*ptInd], 1, 1);
	      book(_h_dijet_cen_groomed[radiusInd][lambdaInd][ptInd], _hepdata_index[11+100*radiusInd+1000*lambdaInd+10000*ptInd], 1, 1);
	      book(_h_dijet_fwd[radiusInd][lambdaInd][ptInd], _hepdata_index[2+100*radiusInd+1000*lambdaInd+10000*ptInd], 1, 1);
	      book(_h_dijet_fwd_groomed[radiusInd][lambdaInd][ptInd], _hepdata_index[12+100*radiusInd+1000*lambdaInd+10000*ptInd], 1, 1);
	    }
	  }
	}
      }
    }

    // Get index of largest bin smaller than value in vector
    // e.g. what you'd need when binning a continuous variable
    uint getBinIndex(float value, const vector<float> & bins) {
      auto itr = std::lower_bound(bins.begin(), bins.end(), value);
      return itr - bins.begin() - 1;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Convert Particles into PseudoJets for clustering
      const VetoedFinalState & fs = apply<VetoedFinalState>(event, "JET_INPUT");
      const Particles & fsParticles = fs.particles();
      vector<PseudoJet> particles;
      particles.reserve(fsParticles.size());
      for (uint iFS=0; iFS<fsParticles.size(); iFS++){
        PseudoJet p = fsParticles[iFS].pseudojet();
        p.set_user_index(fsParticles[iFS].isCharged()); // for later reference to charge
        particles.push_back(p);
      }

      // z jet
      if(_mode==2) {
	for (uint radiusInd=0; radiusInd < _jetRadii.size(); radiusInd++) {
	  float jetRadius = _jetRadii.at(radiusInd);

	  JetDefinition jet_def(antikt_algorithm, jetRadius);
	  vector<PseudoJet> jets = (SelectorPtMin(15))(jet_def(particles));

	  const FinalState& muons = apply<IdentifiedFinalState>(event, "MUONS_NOCUT");
	  if (muons.size() >= 2) {
	    Particle muon1 = muons.particlesByPt()[0];
	    Particle muon2 = muons.particlesByPt()[1];
	    FourMomentum z = muon1.momentum() + muon2.momentum();
	    if (jets.size() > 0) {
	      PseudoJet jet1 = jets[0];
	    }
	  }

	  // Reconstruct Z
	  const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");
	  if (zfinder.bosons().size() < 1) continue;

	  const Particle & z = zfinder.bosons()[0];
	  double zpt = z.pt();

	  // Now do selection criteria
	  bool passZpJ = false;
	  if (jets.size() < 1) continue;
	  PseudoJet jet1 = jets[0];
	  float jet1pt = jet1.pt();
	  float asym = fabs((jet1pt - zpt) / (jet1pt+zpt));
	  float dphi = Rivet::deltaPhi(jet1.phi(), z.phi());
	  passZpJ = ((fabs(jet1.rapidity()) < 1.7) && (zpt > 30) && (asym < 0.3) && (dphi > 2.0));

	  if (!passZpJ) continue;

	  // Now calculate lambda variables and fill hists

	  // Simplify life - ignore this jet if it is below 1st hist pt range
	  // Note that we don't apply it to the original jet pt cut - since
	  // we have phase space where one jet is > 50, and one < 50
	  if (jet1pt < _ptBinsGen[0]) continue;
	  // ignore jet if beyond the last bin
	  if (jet1pt > _ptBinsGen.back()) continue;

	  // Need to use original, ungroomed jet pT to bin
	  uint ptBinInd = getBinIndex(jet1pt, _ptBinsGen);

	  // UNGROOMED VERSION
	  // -------------------------------------------------------------------
	  vector<PseudoJet> chargedParticles;
	  for (uint iC=0; iC<jet1.constituents().size(); iC++)
	    if (jet1.constituents()[iC].user_index())
	      chargedParticles.push_back(jet1.constituents()[iC]);
	  vector<PseudoJet> chargedJets = jet_def(chargedParticles);

	  // Fill hists for each lambda variable
	  for (uint lambdaInd=0; lambdaInd < _lambdaVars.size(); lambdaInd++) {
	    const LambdaVar & thisLambdaVar = _lambdaVars[lambdaInd];
	    Angularity angularity(thisLambdaVar.beta, jetRadius, thisLambdaVar.kappa, thisLambdaVar.constitCut);
	    float val = -1;
	    if (thisLambdaVar.isCharged)
	      val = (chargedJets.size()>0) ? angularity(chargedJets[0]) : -1;
	    else
	      val = angularity(jet1);
	    if (val<0) continue;
	    _h_zpj[radiusInd][lambdaInd][ptBinInd]->fill(val);
	  }

	  // GROOMED VERSION
	  // -------------------------------------------------------------------
	  // Get groomed jet
	  fastjet::contrib::SoftDrop sd(0, 0.1, jetRadius);
	  PseudoJet groomedJet = sd(jet1);
	  PseudoJet groomedJetCharged;
	  if (chargedJets.size()>0)
	    groomedJetCharged= sd(chargedJets[0]);

	  // Fill hists for each lambda variable
	  for (uint lambdaInd=0; lambdaInd < _lambdaVars.size(); lambdaInd++) {
	    const LambdaVar & thisLambdaVar = _lambdaVars[lambdaInd];
	    Angularity angularity(thisLambdaVar.beta, jetRadius, thisLambdaVar.kappa, thisLambdaVar.constitCut);
	    float val = -1;
	    if (thisLambdaVar.isCharged)
	      val = (chargedJets.size()>0) ? angularity(groomedJetCharged) : -1;
	    else
	      val = angularity(groomedJet);
	    if (val<0) continue;
	    _h_zpj_groomed[radiusInd][lambdaInd][ptBinInd]->fill(val);
	  }
	} // end loop over jet radii
      }
      // di jet
      else {
	for (uint radiusInd=0; radiusInd < _jetRadii.size(); radiusInd++) {
	  float jetRadius = _jetRadii.at(radiusInd);

	  JetDefinition jet_def(antikt_algorithm, jetRadius);
	  vector<PseudoJet> jets = (SelectorNHardest(2) * SelectorPtMin(15))(jet_def(particles));

	  bool passDijet = false;
	  if (jets.size() < 2) continue;
	  const auto & jet1 = jets.at(0);
	  const auto & jet2 = jets.at(1);
	  float jet1pt = jet1.pt();
	  float jet2pt = jet2.pt();
	  float asym = (jet1pt - jet2pt) / (jet1pt+jet2pt);
	  float dphi = Rivet::deltaPhi(jet1.phi(), jet2.phi());
	  passDijet = ((fabs(jet1.rapidity()) < 1.7) && (fabs(jet2.rapidity()) < 1.7) && (asym < 0.3) && (dphi > 2.0));

	  if (!passDijet) continue;

	  // Sort by increasing absolute rapidity
	  vector<PseudoJet> dijets = {jet1, jet2};
	  std::sort(dijets.begin(), dijets.end(),
		    [] (const PseudoJet & A, const PseudoJet & B)
		    { return fabs(A.rapidity()) < fabs(B.rapidity()); }
		    );

	  for (uint iJ=0; iJ<dijets.size(); iJ++) {
	    bool isCentral = (iJ == 0);
	    PseudoJet & jetItr = dijets[iJ];

	    // Simplify life - ignore this jet if it is below 1st hist pt range
	    // Note that we don't apply it to the original jet pt cut - since
	    // we have phase space where one jet is > 50, and one < 50
	    if (jetItr.pt() < _ptBinsGen[0]) continue;
	    // ignore jet if beyond the last bin
	    if (jetItr.pt() > _ptBinsGen.back()) continue;

	    // Need to use original, ungroomed jet pT to bin
	    uint ptBinInd = getBinIndex(jetItr.pt(), _ptBinsGen);

	    // UNGROOMED VERSION
	    // -------------------------------------------------------------------
	    vector<PseudoJet> chargedParticles;
	    for (uint iC=0; iC<jetItr.constituents().size(); iC++)
	      if (jetItr.constituents()[iC].user_index())
		chargedParticles.push_back(jetItr.constituents()[iC]);
	    vector<PseudoJet> chargedJets = jet_def(chargedParticles);

	    // Fill hists for each lambda variable
	    for (uint lambdaInd=0; lambdaInd < _lambdaVars.size(); lambdaInd++) {
	      const LambdaVar & thisLambdaVar = _lambdaVars[lambdaInd];
	      Angularity angularity(thisLambdaVar.beta, jetRadius, thisLambdaVar.kappa, thisLambdaVar.constitCut);
	      float val = -1;
	      if (thisLambdaVar.isCharged)
		val = (chargedJets.size()>0) ? angularity(chargedJets[0]) : -1;
	      else
		val = angularity(jetItr);
	      if (val<0) continue;

	      if (isCentral) {
		_h_dijet_cen[radiusInd][lambdaInd][ptBinInd]->fill(val);
	      } else {
		_h_dijet_fwd[radiusInd][lambdaInd][ptBinInd]->fill(val);
	      }
	    }

	    // GROOMED VERSION
	    // -------------------------------------------------------------------
	    // Get groomed jet
	    fastjet::contrib::SoftDrop sd(0, 0.1, jetRadius);
	    PseudoJet groomedJet = sd(jetItr);
	    PseudoJet groomedJetCharged;
	    if (chargedJets.size()>0)
	      groomedJetCharged= sd(chargedJets[0]);

	    // Fill hists for each lambda variable
	    for (uint lambdaInd=0; lambdaInd < _lambdaVars.size(); lambdaInd++) {
	      const LambdaVar & thisLambdaVar = _lambdaVars[lambdaInd];
	      Angularity angularity(thisLambdaVar.beta, jetRadius, thisLambdaVar.kappa, thisLambdaVar.constitCut);
	      float val = -1;
	      if (thisLambdaVar.isCharged)
		val = (chargedJets.size()>0) ? angularity(groomedJetCharged) : -1;
	      else
		val = angularity(groomedJet);
	      if (val<0) continue;

	      if (isCentral) {
		_h_dijet_cen_groomed[radiusInd][lambdaInd][ptBinInd]->fill(val);
	      } else {
		_h_dijet_fwd_groomed[radiusInd][lambdaInd][ptBinInd]->fill(val);
	      }
	    }

	  } // end loop over dijets
	} // end loop over jet radii
      }
    } // end analyze() function


    /// Normalise histograms etc., after the run
    void finalize() {
    } // end of finalize

    /// \class Angularity
    /// definition of angularity
    ///
    class Angularity : public FunctionOfPseudoJet<double>{
    public:
      /// ctor
      Angularity(double alpha, double jet_radius, double kappa=1.0, Selector constitCut=SelectorPtMin(0.)) : _alpha(alpha), _radius(jet_radius), _kappa(kappa), _constitCut(constitCut) {}

      /// computation of the angularity itself
      double result(const PseudoJet &jet) const{
        // get the jet constituents
        vector<PseudoJet> constits = jet.constituents();

        // get the reference axis
        PseudoJet reference_axis = _get_reference_axis(jet);

        // do the actual coputation
        double numerator = 0.0, denominator = 0.0;
        for (const auto &c : constits){
          if (!_constitCut.pass(c)) continue;
          double pt = c.pt();
          // Note: better compute (dist^2)^(alpha/2) to avoid an extra square root
          numerator   += pow(pt, _kappa) * pow(c.squared_distance(reference_axis), 0.5*_alpha);
          denominator += pt;
        }
        if (denominator == 0) return -1;
        // the formula is only correct for the the typical angularities which satisfy either kappa==1 or alpha==0.
        else return numerator/(pow(denominator, _kappa)*pow(_radius, _alpha));
      }

    protected:

      PseudoJet _get_reference_axis(const PseudoJet &jet) const{
        if (_alpha>1) return jet;

        Recluster recluster(JetDefinition(antikt_algorithm, JetDefinition::max_allowable_R, WTA_pt_scheme));
        return recluster(jet);
      }

      double _alpha, _radius, _kappa;
      Selector _constitCut;
    };

    /**
    * Lightweight class to hold info about Lambda variable
    */
    class LambdaVar {

    public:
      LambdaVar(const std::string & name_, float kappa_, float beta_, bool isCharged_, Selector constitCut_):
        name(name_),
        kappa(kappa_),
        beta(beta_),
        isCharged(isCharged_),
        constitCut(constitCut_)
      {}

      std::string name;
      float kappa;
      float beta;
      bool isCharged;
      Selector constitCut;
    };

    // Order matters here
    const vector<float> _jetRadii = {0.4, 0.8};

    // This order is important! index in vector used to create YODA plot name
    // Must match that in extracRivetPlotsDijet.py
    const vector<LambdaVar> _lambdaVars = {
      LambdaVar("jet_multiplicity", 0, 0, false, SelectorPtMin(1.)),
      LambdaVar("jet_pTD", 2, 0, false, SelectorPtMin(0.)),
      LambdaVar("jet_LHA", 1, 0.5, false, SelectorPtMin(0.)),
      LambdaVar("jet_width", 1, 1, false, SelectorPtMin(0.)),
      LambdaVar("jet_thrust", 1, 2, false, SelectorPtMin(0.)),
      LambdaVar("jet_multiplicity_charged", 0, 0, true, SelectorPtMin(1.)),
      LambdaVar("jet_pTD_charged", 2, 0, true, SelectorPtMin(0.)),
      LambdaVar("jet_LHA_charged", 1, 0.5, true, SelectorPtMin(0.)),
      LambdaVar("jet_width_charged", 1, 1, true, SelectorPtMin(0.)),
      LambdaVar("jet_thrust_charged", 1, 2, true, SelectorPtMin(0.)),
    };

    vector<float> _ptBinsGen;

    std::map<unsigned int,unsigned int> _hepdata_index {
{30000,1},{30001,2},{31000,3},{31001,4},{34000,5},{34001,6},{33000,7},{33001,8},{32000,9},{32001,10},{82000,11},{122001,12},{32100,13},{32101,14},{37000,15},{37001,16},{32010,17},{32011,18},{2000,139},{12000,140},{22000,141},{42000,142},{52000,143},{62000,144},{72000,145},{2001,155},{12001,156},{22001,157},{42001,158},{52001,159},{62001,160},{72001,161},{82001,162},{92001,163},{102001,164},{112001,165},{2002,179},{12002,180},{22002,181},{32002,182},{42002,183},{52002,184},{62002,185},{72002,186},{82002,187},{92002,188},{102002,189},{112002,190},{122002,191},{7000,192},{17000,193},{27000,194},{47000,195},{57000,196},{67000,197},{77000,198},{87000,199},{7001,209},{17001,210},{27001,211},{47001,212},{57001,213},{67001,214},{77001,215},{87001,216},{97001,217},{107001,218},{117001,219},{127001,220},{7002,234},{17002,235},{27002,236},{37002,237},{47002,238},{57002,239},{67002,240},{77002,241},{87002,242},{97002,243},{107002,244},{117002,245},{127002,246},{5000,247},{15000,248},{25000,249},{35000,250},{45000,251},{55000,252},{65000,253},{75000,254},{85000,255},{5001,265},{15001,266},{25001,267},{35001,268},{45001,269},{55001,270},{65001,271},{75001,272},{85001,273},{95001,274},{105001,275},{115001,276},{125001,277},{5002,291},{15002,292},{25002,293},{35002,294},{45002,295},{55002,296},{65002,297},{75002,298},{85002,299},{95002,300},{105002,301},{115002,302},{125002,303},{6000,304},{16000,305},{26000,306},{36000,307},{46000,308},{56000,309},{66000,310},{76000,311},{86000,312},{6001,322},{16001,323},{26001,324},{36001,325},{46001,326},{56001,327},{66001,328},{76001,329},{86001,330},{96001,331},{106001,332},{116001,333},{126001,334},{6002,348},{16002,349},{26002,350},{36002,351},{46002,352},{56002,353},{66002,354},{76002,355},{86002,356},{96002,357},{106002,358},{116002,359},{126002,360},{9000,361},{19000,362},{29000,363},{39000,364},{49000,365},{59000,366},{69000,367},{79000,368},{89000,369},{9001,379},{19001,380},{29001,381},{39001,382},{49001,383},{59001,384},{69001,385},{79001,386},{89001,387},{99001,388},{109001,389},{119001,390},{129001,391},{9002,405},{19002,406},{29002,407},{39002,408},{49002,409},{59002,410},{69002,411},{79002,412},{89002,413},{99002,414},{109002,415},{119002,416},{129002,417},{8000,418},{18000,419},{28000,420},{38000,421},{48000,422},{58000,423},{68000,424},{78000,425},{88000,426},{8001,436},{18001,437},{28001,438},{38001,439},{48001,440},{58001,441},{68001,442},{78001,443},{88001,444},{98001,445},{108001,446},{118001,447},{128001,448},{8002,462},{18002,463},{28002,464},{38002,465},{48002,466},{58002,467},{68002,468},{78002,469},{88002,470},{98002,471},{108002,472},{118002,473},{128002,474},{0,475},{10000,476},{20000,477},{40000,478},{50000,479},{60000,480},{70000,481},{80000,482},{1,492},{10001,493},{20001,494},{40001,495},{50001,496},{60001,497},{70001,498},{80001,499},{90001,500},{100001,501},{110001,502},{120001,503},{2,517},{10002,518},{20002,519},{30002,520},{40002,521},{50002,522},{60002,523},{70002,524},{80002,525},{90002,526},{100002,527},{110002,528},{120002,529},{1000,530},{11000,531},{21000,532},{41000,533},{51000,534},{61000,535},{71000,536},{81000,537},{1001,547},{11001,548},{21001,549},{41001,550},{51001,551},{61001,552},{71001,553},{81001,554},{91001,555},{101001,556},{111001,557},{121001,558},{1002,572},{11002,573},{21002,574},{31002,575},{41002,576},{51002,577},{61002,578},{71002,579},{81002,580},{91002,581},{101002,582},{111002,583},{121002,584},{4000,585},{14000,586},{24000,587},{44000,588},{54000,589},{64000,590},{74000,591},{84000,592},{4001,602},{14001,603},{24001,604},{44001,605},{54001,606},{64001,607},{74001,608},{84001,609},{94001,610},{104001,611},{114001,612},{124001,613},{4002,627},{14002,628},{24002,629},{34002,630},{44002,631},{54002,632},{64002,633},{74002,634},{84002,635},{94002,636},{104002,637},{114002,638},{124002,639},{3000,640},{13000,641},{23000,642},{43000,643},{53000,644},{63000,645},{73000,646},{83000,647},{3001,657},{13001,658},{23001,659},{43001,660},{53001,661},{63001,662},{73001,663},{83001,664},{93001,665},{103001,666},{113001,667},{123001,668},{3002,682},{13002,683},{23002,684},{33002,685},{43002,686},{53002,687},{63002,688},{73002,689},{83002,690},{93002,691},{103002,692},{113002,693},{123002,694},{2010,695},{12010,696},{22010,697},{42010,698},{52010,699},{62010,700},{72010,701},{82010,702},{2011,712},{12011,713},{22011,714},{42011,715},{52011,716},{62011,717},{72011,718},{82011,719},{92011,720},{102011,721},{112011,722},{122011,723},{2012,737},{12012,738},{22012,739},{32012,740},{42012,741},{52012,742},{62012,743},{72012,744},{82012,745},{92012,746},{102012,747},{112012,748},{122012,749},{7010,750},{17010,751},{27010,752},{37010,753},{47010,754},{57010,755},{67010,756},{77010,757},{87010,758},{7011,768},{17011,769},{27011,770},{37011,771},{47011,772},{57011,773},{67011,774},{77011,775},{87011,776},{97011,777},{107011,778},{117011,779},{127011,780},{7012,794},{17012,795},{27012,796},{37012,797},{47012,798},{57012,799},{67012,800},{77012,801},{87012,802},{97012,803},{107012,804},{117012,805},{127012,806},{5010,807},{15010,808},{25010,809},{35010,810},{45010,811},{55010,812},{65010,813},{75010,814},{85010,815},{5011,825},{15011,826},{25011,827},{35011,828},{45011,829},{55011,830},{65011,831},{75011,832},{85011,833},{95011,834},{105011,835},{115011,836},{125011,837},{5012,851},{15012,852},{25012,853},{35012,854},{45012,855},{55012,856},{65012,857},{75012,858},{85012,859},{95012,860},{105012,861},{115012,862},{125012,863},{6010,864},{16010,865},{26010,866},{36010,867},{46010,868},{56010,869},{66010,870},{76010,871},{86010,872},{6011,882},{16011,883},{26011,884},{36011,885},{46011,886},{56011,887},{66011,888},{76011,889},{86011,890},{96011,891},{106011,892},{116011,893},{126011,894},{6012,908},{16012,909},{26012,910},{36012,911},{46012,912},{56012,913},{66012,914},{76012,915},{86012,916},{96012,917},{106012,918},{116012,919},{126012,920},{9010,921},{19010,922},{29010,923},{39010,924},{49010,925},{59010,926},{69010,927},{79010,928},{89010,929},{9011,939},{19011,940},{29011,941},{39011,942},{49011,943},{59011,944},{69011,945},{79011,946},{89011,947},{99011,948},{109011,949},{119011,950},{129011,951},{9012,965},{19012,966},{29012,967},{39012,968},{49012,969},{59012,970},{69012,971},{79012,972},{89012,973},{99012,974},{109012,975},{119012,976},{129012,977},{8010,978},{18010,979},{28010,980},{38010,981},{48010,982},{58010,983},{68010,984},{78010,985},{88010,986},{8011,996},{18011,997},{28011,998},{38011,999},{48011,1000},{58011,1001},{68011,1002},{78011,1003},{88011,1004},{98011,1005},{108011,1006},{118011,1007},{128011,1008},{8012,1022},{18012,1023},{28012,1024},{38012,1025},{48012,1026},{58012,1027},{68012,1028},{78012,1029},{88012,1030},{98012,1031},{108012,1032},{118012,1033},{128012,1034},{10,1035},{10010,1036},{20010,1037},{30010,1038},{40010,1039},{50010,1040},{60010,1041},{70010,1042},{80010,1043},{11,1053},{10011,1054},{20011,1055},{30011,1056},{40011,1057},{50011,1058},{60011,1059},{70011,1060},{80011,1061},{90011,1062},{100011,1063},{110011,1064},{120011,1065},{12,1079},{10012,1080},{20012,1081},{30012,1082},{40012,1083},{50012,1084},{60012,1085},{70012,1086},{80012,1087},{90012,1088},{100012,1089},{110012,1090},{120012,1091},{1010,1092},{11010,1093},{21010,1094},{31010,1095},{41010,1096},{51010,1097},{61010,1098},{71010,1099},{81010,1100},{1011,1110},{11011,1111},{21011,1112},{31011,1113},{41011,1114},{51011,1115},{61011,1116},{71011,1117},{81011,1118},{91011,1119},{101011,1120},{111011,1121},{121011,1122},{1012,1136},{11012,1137},{21012,1138},{31012,1139},{41012,1140},{51012,1141},{61012,1142},{71012,1143},{81012,1144},{91012,1145},{101012,1146},{111012,1147},{121012,1148},{4010,1149},{14010,1150},{24010,1151},{34010,1152},{44010,1153},{54010,1154},{64010,1155},{74010,1156},{84010,1157},{4011,1167},{14011,1168},{24011,1169},{34011,1170},{44011,1171},{54011,1172},{64011,1173},{74011,1174},{84011,1175},{94011,1176},{104011,1177},{114011,1178},{124011,1179},{4012,1193},{14012,1194},{24012,1195},{34012,1196},{44012,1197},{54012,1198},{64012,1199},{74012,1200},{84012,1201},{94012,1202},{104012,1203},{114012,1204},{124012,1205},{3010,1206},{13010,1207},{23010,1208},{33010,1209},{43010,1210},{53010,1211},{63010,1212},{73010,1213},{83010,1214},{3011,1224},{13011,1225},{23011,1226},{33011,1227},{43011,1228},{53011,1229},{63011,1230},{73011,1231},{83011,1232},{93011,1233},{103011,1234},{113011,1235},{123011,1236},{3012,1250},{13012,1251},{23012,1252},{33012,1253},{43012,1254},{53012,1255},{63012,1256},{73012,1257},{83012,1258},{93012,1259},{103012,1260},{113012,1261},{123012,1262},{2100,1263},{12100,1264},{22100,1265},{42100,1266},{52100,1267},{62100,1268},{72100,1269},{82100,1270},{2101,1271},{12101,1272},{22101,1273},{42101,1274},{52101,1275},{62101,1276},{72101,1277},{82101,1278},{92101,1279},{102101,1280},{112101,1281},{122101,1282},{2102,1283},{12102,1284},{22102,1285},{32102,1286},{42102,1287},{52102,1288},{62102,1289},{72102,1290},{82102,1291},{92102,1292},{102102,1293},{112102,1294},{122102,1295},{7100,1296},{17100,1297},{27100,1298},{37100,1299},{47100,1300},{57100,1301},{67100,1302},{77100,1303},{87100,1304},{7101,1314},{17101,1315},{27101,1316},{37101,1317},{47101,1318},{57101,1319},{67101,1320},{77101,1321},{87101,1322},{97101,1323},{107101,1324},{117101,1325},{127101,1326},{7102,1327},{17102,1328},{27102,1329},{37102,1330},{47102,1331},{57102,1332},{67102,1333},{77102,1334},{87102,1335},{97102,1336},{107102,1337},{117102,1338},{127102,1339},{5100,1340},{15100,1341},{25100,1342},{35100,1343},{45100,1344},{55100,1345},{65100,1346},{75100,1347},{85100,1348},{5101,1358},{15101,1359},{25101,1360},{35101,1361},{45101,1362},{55101,1363},{65101,1364},{75101,1365},{85101,1366},{95101,1367},{105101,1368},{115101,1369},{125101,1370},{5102,1371},{15102,1372},{25102,1373},{35102,1374},{45102,1375},{55102,1376},{65102,1377},{75102,1378},{85102,1379},{95102,1380},{105102,1381},{115102,1382},{125102,1383},{6100,1384},{16100,1385},{26100,1386},{36100,1387},{46100,1388},{56100,1389},{66100,1390},{76100,1391},{86100,1392},{6101,1402},{16101,1403},{26101,1404},{36101,1405},{46101,1406},{56101,1407},{66101,1408},{76101,1409},{86101,1410},{96101,1411},{106101,1412},{116101,1413},{126101,1414},{6102,1415},{16102,1416},{26102,1417},{36102,1418},{46102,1419},{56102,1420},{66102,1421},{76102,1422},{86102,1423},{96102,1424},{106102,1425},{116102,1426},{126102,1427},{9100,1428},{19100,1429},{29100,1430},{39100,1431},{49100,1432},{59100,1433},{69100,1434},{79100,1435},{89100,1436},{9101,1446},{19101,1447},{29101,1448},{39101,1449},{49101,1450},{59101,1451},{69101,1452},{79101,1453},{89101,1454},{99101,1455},{109101,1456},{119101,1457},{129101,1458},{9102,1459},{19102,1460},{29102,1461},{39102,1462},{49102,1463},{59102,1464},{69102,1465},{79102,1466},{89102,1467},{99102,1468},{109102,1469},{119102,1470},{129102,1471},{8100,1472},{18100,1473},{28100,1474},{38100,1475},{48100,1476},{58100,1477},{68100,1478},{78100,1479},{88100,1480},{8101,1490},{18101,1491},{28101,1492},{38101,1493},{48101,1494},{58101,1495},{68101,1496},{78101,1497},{88101,1498},{98101,1499},{108101,1500},{118101,1501},{128101,1502},{8102,1503},{18102,1504},{28102,1505},{38102,1506},{48102,1507},{58102,1508},{68102,1509},{78102,1510},{88102,1511},{98102,1512},{108102,1513},{118102,1514},{128102,1515},{100,1516},{10100,1517},{20100,1518},{30100,1519},{40100,1520},{50100,1521},{60100,1522},{70100,1523},{80100,1524},{101,1525},{10101,1526},{20101,1527},{30101,1528},{40101,1529},{50101,1530},{60101,1531},{70101,1532},{80101,1533},{90101,1534},{100101,1535},{110101,1536},{120101,1537},{102,1538},{10102,1539},{20102,1540},{30102,1541},{40102,1542},{50102,1543},{60102,1544},{70102,1545},{80102,1546},{90102,1547},{100102,1548},{110102,1549},{120102,1550},{1100,1551},{11100,1552},{21100,1553},{31100,1554},{41100,1555},{51100,1556},{61100,1557},{71100,1558},{81100,1559},{1101,1560},{11101,1561},{21101,1562},{31101,1563},{41101,1564},{51101,1565},{61101,1566},{71101,1567},{81101,1568},{91101,1569},{101101,1570},{111101,1571},{121101,1572},{1102,1573},{11102,1574},{21102,1575},{31102,1576},{41102,1577},{51102,1578},{61102,1579},{71102,1580},{81102,1581},{91102,1582},{101102,1583},{111102,1584},{121102,1585},{4100,1586},{14100,1587},{24100,1588},{34100,1589},{44100,1590},{54100,1591},{64100,1592},{74100,1593},{84100,1594},{4101,1595},{14101,1596},{24101,1597},{34101,1598},{44101,1599},{54101,1600},{64101,1601},{74101,1602},{84101,1603},{94101,1604},{104101,1605},{114101,1606},{124101,1607},{4102,1608},{14102,1609},{24102,1610},{34102,1611},{44102,1612},{54102,1613},{64102,1614},{74102,1615},{84102,1616},{94102,1617},{104102,1618},{114102,1619},{124102,1620},{3100,1621},{13100,1622},{23100,1623},{33100,1624},{43100,1625},{53100,1626},{63100,1627},{73100,1628},{83100,1629},{3101,1630},{13101,1631},{23101,1632},{33101,1633},{43101,1634},{53101,1635},{63101,1636},{73101,1637},{83101,1638},{93101,1639},{103101,1640},{113101,1641},{123101,1642},{3102,1643},{13102,1644},{23102,1645},{33102,1646},{43102,1647},{53102,1648},{63102,1649},{73102,1650},{83102,1651},{93102,1652},{103102,1653},{113102,1654},{123102,1655},{2110,1656},{12110,1657},{22110,1658},{32110,1659},{42110,1660},{52110,1661},{62110,1662},{72110,1663},{82110,1664},{2111,1674},{12111,1675},{22111,1676},{32111,1677},{42111,1678},{52111,1679},{62111,1680},{72111,1681},{82111,1682},{92111,1683},{102111,1684},{112111,1685},{122111,1686},{2112,1687},{12112,1688},{22112,1689},{32112,1690},{42112,1691},{52112,1692},{62112,1693},{72112,1694},{82112,1695},{92112,1696},{102112,1697},{112112,1698},{122112,1699},{7110,1700},{17110,1701},{27110,1702},{37110,1703},{47110,1704},{57110,1705},{67110,1706},{77110,1707},{87110,1708},{7111,1718},{17111,1719},{27111,1720},{37111,1721},{47111,1722},{57111,1723},{67111,1724},{77111,1725},{87111,1726},{97111,1727},{107111,1728},{117111,1729},{127111,1730},{7112,1731},{17112,1732},{27112,1733},{37112,1734},{47112,1735},{57112,1736},{67112,1737},{77112,1738},{87112,1739},{97112,1740},{107112,1741},{117112,1742},{127112,1743},{5110,1744},{15110,1745},{25110,1746},{35110,1747},{45110,1748},{55110,1749},{65110,1750},{75110,1751},{85110,1752},{5111,1762},{15111,1763},{25111,1764},{35111,1765},{45111,1766},{55111,1767},{65111,1768},{75111,1769},{85111,1770},{95111,1771},{105111,1772},{115111,1773},{125111,1774},{5112,1775},{15112,1776},{25112,1777},{35112,1778},{45112,1779},{55112,1780},{65112,1781},{75112,1782},{85112,1783},{95112,1784},{105112,1785},{115112,1786},{125112,1787},{6110,1788},{16110,1789},{26110,1790},{36110,1791},{46110,1792},{56110,1793},{66110,1794},{76110,1795},{86110,1796},{6111,1806},{16111,1807},{26111,1808},{36111,1809},{46111,1810},{56111,1811},{66111,1812},{76111,1813},{86111,1814},{96111,1815},{106111,1816},{116111,1817},{126111,1818},{6112,1819},{16112,1820},{26112,1821},{36112,1822},{46112,1823},{56112,1824},{66112,1825},{76112,1826},{86112,1827},{96112,1828},{106112,1829},{116112,1830},{126112,1831},{9110,1832},{19110,1833},{29110,1834},{39110,1835},{49110,1836},{59110,1837},{69110,1838},{79110,1839},{89110,1840},{9111,1850},{19111,1851},{29111,1852},{39111,1853},{49111,1854},{59111,1855},{69111,1856},{79111,1857},{89111,1858},{99111,1859},{109111,1860},{119111,1861},{129111,1862},{9112,1863},{19112,1864},{29112,1865},{39112,1866},{49112,1867},{59112,1868},{69112,1869},{79112,1870},{89112,1871},{99112,1872},{109112,1873},{119112,1874},{129112,1875},{8110,1876},{18110,1877},{28110,1878},{38110,1879},{48110,1880},{58110,1881},{68110,1882},{78110,1883},{88110,1884},{8111,1894},{18111,1895},{28111,1896},{38111,1897},{48111,1898},{58111,1899},{68111,1900},{78111,1901},{88111,1902},{98111,1903},{108111,1904},{118111,1905},{128111,1906},{8112,1907},{18112,1908},{28112,1909},{38112,1910},{48112,1911},{58112,1912},{68112,1913},{78112,1914},{88112,1915},{98112,1916},{108112,1917},{118112,1918},{128112,1919},{110,1920},{10110,1921},{20110,1922},{30110,1923},{40110,1924},{50110,1925},{60110,1926},{70110,1927},{80110,1928},{111,1938},{10111,1939},{20111,1940},{30111,1941},{40111,1942},{50111,1943},{60111,1944},{70111,1945},{80111,1946},{90111,1947},{100111,1948},{110111,1949},{120111,1950},{112,1951},{10112,1952},{20112,1953},{30112,1954},{40112,1955},{50112,1956},{60112,1957},{70112,1958},{80112,1959},{90112,1960},{100112,1961},{110112,1962},{120112,1963},{1110,1964},{11110,1965},{21110,1966},{31110,1967},{41110,1968},{51110,1969},{61110,1970},{71110,1971},{81110,1972},{1111,1982},{11111,1983},{21111,1984},{31111,1985},{41111,1986},{51111,1987},{61111,1988},{71111,1989},{81111,1990},{91111,1991},{101111,1992},{111111,1993},{121111,1994},{1112,1995},{11112,1996},{21112,1997},{31112,1998},{41112,1999},{51112,2000},{61112,2001},{71112,2002},{81112,2003},{91112,2004},{101112,2005},{111112,2006},{121112,2007},{4110,2008},{14110,2009},{24110,2010},{34110,2011},{44110,2012},{54110,2013},{64110,2014},{74110,2015},{84110,2016},{4111,2026},{14111,2027},{24111,2028},{34111,2029},{44111,2030},{54111,2031},{64111,2032},{74111,2033},{84111,2034},{94111,2035},{104111,2036},{114111,2037},{124111,2038},{4112,2039},{14112,2040},{24112,2041},{34112,2042},{44112,2043},{54112,2044},{64112,2045},{74112,2046},{84112,2047},{94112,2048},{104112,2049},{114112,2050},{124112,2051},{3110,2052},{13110,2053},{23110,2054},{33110,2055},{43110,2056},{53110,2057},{63110,2058},{73110,2059},{83110,2060},{3111,2070},{13111,2071},{23111,2072},{33111,2073},{43111,2074},{53111,2075},{63111,2076},{73111,2077},{83111,2078},{93111,2079},{103111,2080},{113111,2081},{123111,2082},{3112,2083},{13112,2084},{23112,2085},{33112,2086},{43112,2087},{53112,2088},{63112,2089},{73112,2090},{83112,2091},{93112,2092},{103112,2093},{113112,2094},{123112,2095},};

    // mode for the analysis
    unsigned int _mode;
    // 3D vector: [jet radius][lambda variable][pt bin]
    // since each pt bin has its own normalised distribution
    vector<vector<vector<Histo1DPtr> > > _h_dijet_cen,
                                         _h_dijet_cen_groomed,
                                         _h_dijet_fwd,
                                         _h_dijet_fwd_groomed;
    // 3D vector: [jet radius][lambda variable][pt bin]
    // since each pt bin has its own normalised distribution
    vector<vector<vector<Histo1DPtr> > > _h_zpj, _h_zpj_groomed;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2021_I1920187);

}