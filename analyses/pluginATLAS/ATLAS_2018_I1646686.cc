// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Math/LorentzTrans.hh"

namespace Rivet {

  class ATLAS_2018_I1646686 : public Analysis {

    public:
      /// Constructor
      /// @brief all-hadronic ttbar at 13 TeV
      DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2018_I1646686);

      /// Book histograms and initialise projections before the run
      void init() {

        //histogram booking
        _h["inclusive"] = bookHisto1D(2,1,1);
        bookHistograms("t_y",  1, true);
        bookHistograms("t_pt", 2, true);
        bookHistograms("t1_pt",         3);
        bookHistograms("t1_y",          4);
        bookHistograms("t2_pt",         5);
        bookHistograms("t2_y",          6);
        bookHistograms("tt_m",          7);
        bookHistograms("tt_pt",         8);
        bookHistograms("tt_y",          9);
        bookHistograms("tt_chi",       10);
        bookHistograms("tt_yboost",    11);
        bookHistograms("tt_pout",      12);
        bookHistograms("tt_dPhi",      13);
        bookHistograms("tt_Ht",        14);
        bookHistograms("tt_cosThStar", 15);

        //projections
        Cut dressed_lep = (Cuts::abseta < 2.5) && (Cuts::pT >= 25*GeV);
        Cut eta_full = (Cuts::abseta < 5.0);

        // All final state particles
        FinalState fs(eta_full);

        // Get photons to dress leptons
        IdentifiedFinalState photons(fs);
        photons.acceptIdPair(PID::PHOTON);

        // Projection to find the electrons
        IdentifiedFinalState el_id(fs);
        el_id.acceptIdPair(PID::ELECTRON);
        PromptFinalState electrons(el_id);
        electrons.acceptTauDecays(true);
        DressedLeptons dressedelectrons(photons, electrons, 0.1, dressed_lep);
        addProjection(dressedelectrons, "elecs");
        DressedLeptons ewdressedelectrons(photons, electrons, 0.1, eta_full);

        // Projection to find the muons
        IdentifiedFinalState mu_id(fs);
        mu_id.acceptIdPair(PID::MUON);
        PromptFinalState muons(mu_id);
        muons.acceptTauDecays(true);
        DressedLeptons dressedmuons(photons, muons, 0.1, dressed_lep);
        addProjection(dressedmuons, "muons");
        DressedLeptons ewdressedmuons(photons, muons, 0.1, eta_full);

        // Projection to find neutrinos for vetoing in jets 
        IdentifiedFinalState nu_id(fs);
        nu_id.acceptNeutrinos();
        PromptFinalState neutrinos(nu_id);
        neutrinos.acceptTauDecays(true);

        // Jet clustering.
        VetoedFinalState lvfs(fs);
        lvfs.addVetoOnThisFinalState(mu_id);
        lvfs.addVetoOnThisFinalState(nu_id);

        VetoedFinalState vfs;
        vfs.addVetoOnThisFinalState(ewdressedelectrons);
        vfs.addVetoOnThisFinalState(ewdressedmuons);
        vfs.addVetoOnThisFinalState(neutrinos);
  
        FastJets jets(vfs, FastJets::ANTIKT, 0.4);
        jets.useInvisibles();
        addProjection(jets, "jets");

        FastJets ljets(lvfs, FastJets::ANTIKT, 1.0);
        ljets.useInvisibles(); 
        addProjection(ljets, "ljets" );

        PartonicTops partonTops;
        addProjection(partonTops, "partonicTops");
      }

      //////////////////////////

      void analyze(const Event& event) {

        const double weight = event.weight();

        const Particles partonicTops = applyProjection<PartonicTops>( event, "partonicTops").particlesByPt();

        // Parton-level top quarks
        FourMomentum top;
        FourMomentum tbar;
        bool foundT = false;
        bool foundTBar = false;
        for (const Particle& ptop : partonicTops) {
          const int pid = ptop.pdgId();
          if (pid == 6) {
            top = ptop.momentum();
            foundT = true;
          } 
          else if (pid == -6) { 
            tbar = ptop.momentum();
            foundTBar = true;
          }
        } //loop over partonic tops

        FourMomentum t1_parton;
        FourMomentum t2_parton;
        FourMomentum ttbar_parton;

        if ( foundT && foundTBar ) {
          t1_parton = top.pT2() > tbar.pT2() ? top : tbar;
          t2_parton = top.pT2() > tbar.pT2() ? tbar : top;
          ttbar_parton = t1_parton + t2_parton;

          if ( t1_parton.pT() > 500*GeV && t2_parton.pT() > 350*GeV) {

            double chi_parton = CalcChi(t1_parton, t2_parton);
            double cosThetaStar_parton = abs(CalcCosThetaStar(t1_parton, t2_parton));
            double pout_parton = abs(CalcPout(t1_parton, t2_parton));
            double dPhi_parton = deltaPhi(t1_parton, t2_parton);

            int randomChoice = rand() % 2;
            FourMomentum& randomTopParton = (randomChoice == 0) ? t1_parton : t2_parton;
            fillParton("t_pt", randomTopParton.pT()/GeV, weight);
            fillParton("t_y",  randomTopParton.absrap(), weight);

            fillParton("t1_pt", t1_parton.pT()/GeV, weight);
            fillParton("t1_y",  t1_parton.absrap(), weight);
            fillParton("t2_pt", t2_parton.pT()/GeV, weight);
            fillParton("t2_y",  t2_parton.absrap(), weight);

            fillParton("tt_m",  ttbar_parton.mass()/GeV, weight);
            fillParton("tt_pt", ttbar_parton.pT()/GeV, weight);
            fillParton("tt_Ht", (t1_parton.pT() + t2_parton.pT())/GeV, weight);
            fillParton("tt_y",  ttbar_parton.absrap(), weight);

            fillParton("tt_yboost", 0.5 * abs(t1_parton.rapidity() + t2_parton.rapidity()), weight);
            fillParton("tt_chi", chi_parton, weight);
            fillParton("tt_cosThStar", cosThetaStar_parton, weight);
            fillParton("tt_pout", pout_parton/GeV, weight);
            fillParton("tt_dPhi", dPhi_parton, weight);
          } 
        } //found partonic tops

        const vector<DressedLepton> dressedElectrons = applyProjection<DressedLeptons>(event, "elecs").dressedLeptons();
        const vector<DressedLepton> dressedMuons     = applyProjection<DressedLeptons>(event, "muons").dressedLeptons();

        const Jets& all_jets  = applyProjection<FastJets>( event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);

        const FastJets& ljets_fj = applyProjection<FastJets>( event, "ljets");
        const Jets all_ljets    = ljets_fj.jetsByPt();

        PseudoJets trimmed_fat_pjets;
        Jets trimmedJets;
        fastjet::Filter trimmer(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.05));
        for (const Jet& jet : all_ljets)  trimmedJets += ljets_fj.trimJet(jet, trimmer);
        trimmedJets = Rivet::sortByPt(trimmedJets);  

        if ( dressedElectrons.size() )   vetoEvent;
        else if ( dressedMuons.size() )  vetoEvent;

        // Check largeR jets
        Jets ljets;
        std::vector<bool> b_tagged; 
        for (const Jet& jet : trimmedJets) {

          if (jet.pT() < 250 * GeV)  continue; 
          if (jet.pT() > 3000 * GeV) continue; 
          if (jet.mass() > jet.pT()) continue; 
          if (jet.abseta() > 2.0 )   continue;

          ljets += jet;
          bool btagged = !jet.bTags().empty();
          b_tagged.push_back(btagged);
        } // loop over ljets

        if (all_jets.size() < 2)  vetoEvent;
        if (ljets.size() < 2)     vetoEvent;

        //identify top and anti top, compute some event variables
        FourMomentum ttbar = ljets[0].momentum() + ljets[1].momentum();                                
        FourMomentum t1 = ljets[0].momentum(); 
        FourMomentum t2 = ljets[1].momentum(); 

        double chi = CalcChi(t1, t2);
        double cosThetaStar = abs(CalcCosThetaStar(t1, t2));
        double pout = abs(CalcPout(t1, t2));
        double dPhi = deltaPhi(t1, t2);

        if ( t2.pT() < 350*GeV)  vetoEvent;
        if ( t1.pT() < 500*GeV)  vetoEvent;

        //b-tagging for particle done on large-R jets
        if (!(b_tagged[0] && b_tagged[1]))  vetoEvent; 

        //continues with signal region cuts
        if ( abs(t1.mass() - 172.5 * GeV) > 50*GeV )  vetoEvent;
        if ( abs(t2.mass() - 172.5 * GeV) > 50*GeV )  vetoEvent;

        _h["inclusive"]->fill(0, weight);  

        fillHistograms("t1_pt", t1.pT()/GeV, weight);
        fillHistograms("t1_y",  t1.absrap(), weight);
        fillHistograms("t2_pt", t2.pT()/GeV, weight);
        fillHistograms("t2_y",  t2.absrap(), weight);

        fillHistograms("tt_m",  ttbar.mass()/TeV, weight);
        fillHistograms("tt_pt", ttbar.pT()/GeV, weight);
        fillHistograms("tt_Ht", (t1.pT() + t2.pT())/GeV, weight);
        fillHistograms("tt_y",  ttbar.absrap(), weight);

        fillHistograms("tt_yboost", 0.5 * abs(t1.rapidity() + t2.rapidity()), weight);
        fillHistograms("tt_chi", chi, weight);
        fillHistograms("tt_cosThStar", cosThetaStar, weight);
        fillHistograms("tt_pout", pout/GeV, weight);
        fillHistograms("tt_dPhi", dPhi, weight);

      } //analyze

      //////////////////////////////////////////////////

      void finalize() {
        // Normalize histograms
        const double sf = crossSection() / sumOfWeights();
        for (auto &hist : _h) {
          scale(hist.second, sf);
          if (hist.first.find("_norm") != string::npos)  hist.second->normalize(1.0, false);
        } 
      } // finalize

      double CalcChi(const FourMomentum& t1, const FourMomentum& t2) {   
        double ystar = 0.5 * (t1.rapidity()-t2.rapidity());
        double chi = exp( 2 * abs(ystar)); 
        return chi;
      } //CalcChi

      double CalcCosThetaStar(const FourMomentum& t1, const FourMomentum& t2) {
        FourMomentum ttbar = t1 + t2;
        LorentzTransform centreOfMassTrans;
        ttbar.setX(0);
        ttbar.setY(0);
        centreOfMassTrans.setBetaVec( -ttbar.boostVector() );
        FourMomentum t1_star = centreOfMassTrans.transform(t1);
        double cosThetaStar = t1_star.pz()/t1_star.p3().mod(); 
        return cosThetaStar;
      } //CalcCosThetaStar

      double CalcPout(const FourMomentum& t1, const FourMomentum& t2) { 
        Vector3 t1V = t1.p3();
        Vector3 t2V = t2.p3();
        Vector3 zUnit = Vector3(0., 0., 1.);
        Vector3 vPerp = zUnit.cross(t1V);

        double pout = vPerp.dot(t2V)/vPerp.mod(); 
        return pout;
      } //CalcPout

    private:

      map<string, Histo1DPtr> _h;

      //some functions for booking, filling and scaling the histograms
      void fillHistograms(std::string name, double value, double weight) { 
        _h[name]->fill(value, weight);
        _h[name + "_norm"]->fill(value, weight);
      } 

      void fillParton(std::string name, double value, double weight) { 
        _h[name + "_parton"]->fill(value, weight);
        _h[name + "_parton_norm"]->fill(value, weight);
      } 

      void bookHistograms(const std::string name, unsigned int index, bool onlyParton = false) {
        if (!onlyParton) {
          _h[name] = bookHisto1D(index, 1, 1 ); 
          _h[name + "_norm"] = bookHisto1D(index + 13, 1, 1 ); 
        }
        _h[name + "_parton"] = bookHisto1D(index + 82, 1, 1 ); 
        _h[name + "_parton_norm"] = bookHisto1D(index + 97, 1, 1 ); 
      } 

  }; // end class ATLAS_2018_I1646686

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2018_I1646686);
}
