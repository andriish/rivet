#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/ParticleName.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {

  class CMS_2016_I1473674 : public Analysis {
  public:

    // Minimal constructor
    CMS_2016_I1473674() : Analysis("CMS_2016_I1473674") {
    }

    // Set up projections and book histograms
    void init() {
      // Complete final state
      FinalState fs(-MAXDOUBLE, MAXDOUBLE, 0*GeV);

      // Projection for taus
      TauFinder taus(TauFinder::ANY);
      addProjection(taus, "Tau");
      IdentifiedFinalState nu_taus(fs);
      nu_taus.acceptIdPair(PID::NU_TAU);
      addProjection(nu_taus, "NuTau");

      // Projection for electrons and muons
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);

      IdentifiedFinalState el_id(fs);
      el_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState electrons(el_id);
      addProjection(electrons, "Electrons");
      DressedLeptons dressed_electrons(photons, electrons, 0.1, Cuts::open(), true, false);
      addProjection(dressed_electrons, "DressedElectrons");

      IdentifiedFinalState mu_id(fs);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState muons(mu_id);
      addProjection(muons, "Muons");
      DressedLeptons dressed_muons(photons, muons, 0.1, Cuts::open(), true, false);
      addProjection(dressed_muons, "DressedMuons");

      // Projection for jets
      VetoedFinalState fs_jets(FinalState(-MAXDOUBLE, MAXDOUBLE, 0*GeV));
      fs_jets.addVetoOnThisFinalState(dressed_muons);
      addProjection(FastJets(fs_jets, FastJets::ANTIKT, 0.5), "Jets");

      // Projections for MET
      addProjection(MissingMomentum(), "MET");

      // Booking of histograms
      _hist_met = bookHisto1D(5, 1, 1);
      _hist_ht  = bookHisto1D(6, 1, 1);
      _hist_st  = bookHisto1D(7, 1, 1);
      _hist_wpt = bookHisto1D(8, 1, 1);
    }


    // per event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // select ttbar -> lepton+jets
      const DressedLeptons& dressed_electrons = applyProjection<DressedLeptons>(event, "DressedElectrons");
      const DressedLeptons& dressed_muons = applyProjection<DressedLeptons>(event, "DressedMuons");
      if (dressed_electrons.dressedLeptons().size() +
          dressed_muons.dressedLeptons().size() != 1) {
        vetoEvent;
      }

      FourMomentum lepton;
      if (dressed_electrons.dressedLeptons().size() == 1) {
        lepton = dressed_electrons.dressedLeptons()[0].momentum();
      } else {
        lepton = dressed_muons.dressedLeptons()[0].momentum();
      }

      // veto if lepton is tau
      const TauFinder& taus = applyProjection<TauFinder>(event, "Tau");
      const IdentifiedFinalState nu_taus = applyProjection<IdentifiedFinalState>(event, "NuTau");
      foreach (const Particle& tau, taus.taus()) {
        foreach (const Particle& nu, nu_taus.particles()) {
          if (tau.pid() * nu.pid() < 0)
            continue;

          const FourMomentum w_candidate = tau.momentum() + nu.momentum();
          if (abs(w_candidate.mass() - 80.4) > 5.)
            vetoEvent;
        }
      }


      // MET
      const MissingMomentum& met = applyProjection<MissingMomentum>(event, "MET");
      _hist_met->fill(met.visibleMomentum().pT()/GeV, weight);

      // HT and ST
      const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
      const Jets jets = jetpro.jetsByPt(20*GeV);

      double ht = 0.0;
      foreach (const Jet& j, jets) {
        if (deltaR(j.momentum(), lepton) > 0.3) {
          ht += j.pT();
        }
      }

      double st = ht + lepton.pT() + met.visibleMomentum().pT();
      _hist_ht->fill(ht/GeV, weight);
      _hist_st->fill(st/GeV, weight);

      // WPT
      FourMomentum w = lepton - met.visibleMomentum();
      _hist_wpt->fill(w.pT()/GeV, weight);
    }

    // scale by 1 over weight
    void finalize() {
      normalize(_hist_met);
      normalize(_hist_ht);
      normalize(_hist_st);
      normalize(_hist_wpt);
    }

  private:
    Histo1DPtr _hist_met, _hist_ht, _hist_st, _hist_wpt;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_I1473674);
}
