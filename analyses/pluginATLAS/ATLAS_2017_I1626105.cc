#include "Rivet/Analysis.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// @brief lepton differential ttbar analysis at 8 TeV
  class ATLAS_2017_I1626105 : public Analysis {
  public:
    
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2017_I1626105);
    
    void init() {

      Cut eta_full = Cuts::abseta < 5.0 && Cuts::pT > 1.0*MeV;

      // All final state particles
      const FinalState fs;

      // Get photons to dress leptons
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);

      // Projection to find the electrons
      PromptFinalState prompt_el(Cuts::abspid == PID::ELECTRON, true);
      DressedLeptons elecs(photons, prompt_el, 0.1, (Cuts::abseta < 2.5) && (Cuts::pT > 25*GeV));
      DressedLeptons veto_elecs(photons, prompt_el, 0.1, eta_full, false);
      declare(elecs, "elecs");

      // Projection to find the muons
      PromptFinalState prompt_mu(Cuts::abspid == PID::MUON, true);
      DressedLeptons muons(photons, prompt_mu, 0.1, (Cuts::abseta < 2.5) && (Cuts::pT > 25*GeV));
      DressedLeptons veto_muons(photons, prompt_mu, 0.1, eta_full, false);
      declare(muons, "muons");

      // Jet clustering.
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(veto_elecs);
      vfs.addVetoOnThisFinalState(veto_muons);
      FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles();
      declare(jets, "jets");

      // Book histograms
      bookHistos("lep_pt",       1);
      bookHistos("lep_eta",      3);
      bookHistos("dilep_pt",     5);
      bookHistos("dilep_mass",   7);
      bookHistos("dilep_rap",    9);
      bookHistos("dilep_dphi",  11);
      bookHistos("dilep_sumpt", 13);
      bookHistos("dilep_sumE",  15);
    }

    void analyze(const Event& event) {
      vector<DressedLepton> elecs = sortByPt(apply<DressedLeptons>(event, "elecs").dressedLeptons());
      vector<DressedLepton> muons = sortByPt(apply<DressedLeptons>(event, "muons").dressedLeptons());      
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
     
      // Check overlap of jets/leptons.
      for (const Jet& jet : jets) {
        for (const DressedLepton& el : elecs) {
          if (deltaR(jet, el) < 0.4) delete el;
        }
        for (const DressedLepton& mu : muons) {
          if (deltaR(jet, mu) < 0.4) delete mu;
        }
      }
      if (elecs.empty() || muons.empty())  vetoEvent;
      
      if (elecs[0].charge() == muons[0].charge())  vetoEvent;  
      
      FourMomentum el = elecs[0].momentum();
      FourMomentum mu = muons[0].momentum();
      FourMomentum ll = elecs[0].momentum() + muons[0].momentum();
                  
      // Fill histograms
      const double weight = event.weight();
      fillHistos("lep_pt",      el.pT()/GeV,             weight);
      fillHistos("lep_pt",      mu.pT()/GeV,             weight);
      fillHistos("lep_eta",     el.abseta(),             weight);
      fillHistos("lep_eta",     mu.abseta(),             weight);
      fillHistos("dilep_pt",    ll.pT()/GeV,             weight);
      fillHistos("dilep_mass",  ll.mass()/GeV,           weight);
      fillHistos("dilep_rap",   ll.absrap(),             weight);
      fillHistos("dilep_dphi",  deltaPhi(el, mu),        weight);
      fillHistos("dilep_sumpt", (el.pT() + mu.pT())/GeV, weight);
      fillHistos("dilep_sumE",  (el.E() + mu.E())/GeV,   weight);
    }

    void finalize() {
      // Normalize to cross-section
      const double sf = crossSection()/femtobarn/sumOfWeights();
      for (auto& hist : _h) {
        const double norm = 1.0 / hist.second->integral();
        // add overflow to last bin 
        double overflow = hist.second->overflow().effNumEntries();
        hist.second->fillBin(hist.second->numBins() - 1, overflow);
        // histogram normalisation
        if (hist.first.find("norm") != string::npos)  scale(hist.second, norm);
        else  scale(hist.second, sf);
      }

    }

  private:

    /// @name Histogram helper functions
    //@{
    void bookHistos(const std::string name, unsigned int index) {
      _h[name] = bookHisto1D(index, 1, 1);
      _h["norm_" + name] = bookHisto1D(index + 1, 1, 1);
    }

    void fillHistos(const std::string name, double value, double weight) {
      _h[name]->fill(value, weight);
      _h["norm_" + name]->fill(value, weight);
    }

    map<string, Histo1DPtr> _h;
    //@}

  };

  // Declare the class as a hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1626105);
}
