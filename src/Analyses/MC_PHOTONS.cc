// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class MC_PHOTONS : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    MC_PHOTONS()
      : Analysis("MC_PHOTONS")
    {    }

    //@}


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      IdentifiedFinalState leptons(-5.0, 5.0, 10*GeV);
      IdentifiedFinalState photons(-5.0, 5.0);
      leptons.acceptChLeptons();
      photons.acceptId(PHOTON);
      addProjection(leptons, "lFS");
      addProjection(photons, "gammaFS");

      _h_Ptgamma = bookHistogram1D("Ptgamma", logspace(50, 0.01, 30));
      _h_Egamma = bookHistogram1D("Egamma", logspace(50, 0.01, 200));
      _h_sumPtgamma = bookHistogram1D("sumPtgamma", logspace(50, 0.1, 100));
      _h_sumEgamma = bookHistogram1D("sumEgamma", logspace(50, 0.1, 500));
      _h_DelR = bookHistogram1D("DeltaR", 50, 0, 2);
      _h_DelR_weighted = bookHistogram1D("DeltaR_weighted", 50, 0, 2);
      _h_DelR_R = bookHistogram1D("DeltaR_R", 50, 0, 2);
      _h_DelR_R_weighted = bookHistogram1D("DeltaR_R_weighted", 50, 0, 2);
      _p_DelR_vs_pTl = bookProfile1D("DeltaR_vs_pTlep", logspace(50, 10, 50*log10(sqrtS())));
      _p_DelR_weighted_vs_pTl = bookProfile1D("DeltaR_weighted_vs_pTlep", logspace(50, 10, 50*log10(sqrtS())));
      _p_sumPtgamma_vs_pTl = bookProfile1D("sumPtGamma_vs_pTlep", logspace(50, 10, 50*log10(sqrtS())));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      /// Get photons and leptons
      const ParticleVector& photons = applyProjection<FinalState>(event, "gammaFS").particles();
      MSG_DEBUG("Photon multiplicity = " << photons.size());
      const ParticleVector& leptons = applyProjection<FinalState>(event, "lFS").particles();
      MSG_DEBUG("Photon multiplicity = " << leptons.size());

      // Initialise a map of sumPtgamma for each lepton
      map<size_t, double> sumpT_per_lep;
      for (size_t il = 0; il < leptons.size(); ++il) sumpT_per_lep[il] = 0;

      // Calculate photon energies and transverse momenta
      double sumPtgamma(0), sumEgamma(0);
      foreach (const Particle& p, photons) {
        // Individual and summed pTs and energies
        double pTgamma = p.momentum().pT()/GeV;
        double Egamma = p.momentum().E()/GeV;
        _h_Ptgamma->fill(pTgamma, weight);
        _h_Egamma->fill(Egamma, weight);
        sumPtgamma += pTgamma;
        sumEgamma += Egamma;

        // Calculate delta R with respect to the nearest lepton
        int ilep = -1;
        double delR = 10000;
        for (size_t il = 0; il < leptons.size(); ++il) {
          const double tmpdelR = deltaR(leptons[il].momentum(), p.momentum());
          if (tmpdelR < delR) {
            ilep = il;
            delR = tmpdelR;
          }
        }
        if (ilep != -1) {
          _h_DelR->fill(delR, weight);
          _h_DelR_weighted->fill(delR*pTgamma/GeV, weight);
          _h_DelR_R->fill(delR, weight/(delR+1e-5));
          _h_DelR_R_weighted->fill(delR*pTgamma/GeV, weight/(delR+1e-5));
          _p_DelR_vs_pTl->fill(leptons[ilep].momentum().pT()/GeV, delR, weight);
          _p_DelR_weighted_vs_pTl->fill(leptons[ilep].momentum().pT()/GeV, delR*pTgamma, weight);
          sumpT_per_lep[ilep] += pTgamma;
        }
      }

      // Histogram whole-event photon HT/energy
      _h_sumPtgamma->fill(sumPtgamma/GeV, weight);
      _h_sumEgamma->fill(sumEgamma/GeV, weight);

      // Histogram per-lepton sum(pT)
      for (size_t il = 0; il < leptons.size(); ++il) {
        _p_sumPtgamma_vs_pTl->fill(leptons[il].momentum().pT()/GeV, sumpT_per_lep[il]/GeV, weight);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Ptgamma);
      normalize(_h_Egamma);
      normalize(_h_sumPtgamma);
      normalize(_h_sumEgamma);
      normalize(_h_DelR);
      normalize(_h_DelR_weighted);
      normalize(_h_DelR_R);
      normalize(_h_DelR_R_weighted);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D *_h_Ptgamma, *_h_Egamma;
    AIDA::IHistogram1D *_h_sumPtgamma, *_h_sumEgamma;
    AIDA::IHistogram1D *_h_DelR, *_h_DelR_weighted;
    AIDA::IHistogram1D *_h_DelR_R, *_h_DelR_R_weighted;
    AIDA::IProfile1D *_p_DelR_vs_pTl, *_p_DelR_weighted_vs_pTl;
    AIDA::IProfile1D *_p_sumPtgamma_vs_pTl;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_PHOTONS);


}
