// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Particle.fhh"

namespace Rivet {


  class ATLAS_2012_I1204784 : public Analysis {
    public:

      /// Constructor
      ATLAS_2012_I1204784()
        : Analysis("ATLAS_2012_I1204784")
      {      }


    public:

      /// Book histograms and initialise projections before the run
      void init() {

        ZFinder zfinder_dressed_el(-2.4, 2.4, 20, ELECTRON, 66.0*GeV, 116.0*GeV, 0.1, true, false);
        addProjection(zfinder_dressed_el, "ZFinder_dressed_el");
        ZFinder zfinder_bare_el(-2.4, 2.4, 20, ELECTRON, 66.0*GeV, 116.0*GeV, 0.0, true, false);
        addProjection(zfinder_bare_el, "ZFinder_bare_el");
        ZFinder zfinder_dressed_mu(-2.4, 2.4, 20, MUON, 66.0*GeV, 116.0*GeV, 0.1, true, false);
        addProjection(zfinder_dressed_mu, "ZFinder_dressed_mu");
        ZFinder zfinder_bare_mu(-2.4, 2.4, 20, MUON, 66.0*GeV, 116.0*GeV, 0.0, true, false);
        addProjection(zfinder_bare_mu, "ZFinder_bare_mu");

        // Book histograms
        // Single-differential plots
        _hist_zphistar_el_bare = bookHistogram1D(1, 1, 1);
        _hist_zphistar_mu_bare = bookHistogram1D(1, 1, 2);
        _hist_zphistar_el_dressed = bookHistogram1D(2, 1, 1);
        _hist_zphistar_mu_dressed = bookHistogram1D(2, 1, 2);

        // Double-differential plots
        _h_phistar_el_bare.addHistogram(0.0, 0.8, bookHistogram1D(3, 1, 1));
        _h_phistar_el_bare.addHistogram(0.8, 1.6, bookHistogram1D(3, 1, 2));
        _h_phistar_el_bare.addHistogram(1.6, 10.0, bookHistogram1D(3, 1, 3));

        _h_phistar_el_dressed.addHistogram(0.0, 0.8, bookHistogram1D(3, 2, 1));
        _h_phistar_el_dressed.addHistogram(0.8, 1.6, bookHistogram1D(3, 2, 2));
        _h_phistar_el_dressed.addHistogram(1.6, 10.0, bookHistogram1D(3, 2, 3));

        _h_phistar_mu_bare.addHistogram(0.0, 0.8, bookHistogram1D(4, 1, 1));
        _h_phistar_mu_bare.addHistogram(0.8, 1.6, bookHistogram1D(4, 1, 2));
        _h_phistar_mu_bare.addHistogram(1.6, 10.0, bookHistogram1D(4, 1, 3));

        _h_phistar_mu_dressed.addHistogram(0.0, 0.8, bookHistogram1D(4, 2, 1));
        _h_phistar_mu_dressed.addHistogram(0.8, 1.6, bookHistogram1D(4, 2, 2));
        _h_phistar_mu_dressed.addHistogram(1.6, 10.0, bookHistogram1D(4, 2, 3));
      }


      /// Perform the per-event analysis
      void analyze(const Event& event) {
        const double weight = event.weight();

        const ZFinder& zfinder_dressed_el = applyProjection<ZFinder>(event, "ZFinder_dressed_el");
        const ZFinder& zfinder_bare_el = applyProjection<ZFinder>(event, "ZFinder_bare_el");
        const ZFinder& zfinder_dressed_mu = applyProjection<ZFinder>(event, "ZFinder_dressed_mu");
        const ZFinder& zfinder_bare_mu = applyProjection<ZFinder>(event, "ZFinder_bare_mu");

        fillPlots(zfinder_dressed_el, _hist_zphistar_el_dressed, _h_phistar_el_dressed, weight);
        fillPlots(zfinder_bare_el, _hist_zphistar_el_bare, _h_phistar_el_bare, weight);
        fillPlots(zfinder_dressed_mu, _hist_zphistar_mu_dressed, _h_phistar_mu_dressed, weight);
        fillPlots(zfinder_bare_mu, _hist_zphistar_mu_bare, _h_phistar_mu_bare, weight);
      }


      void fillPlots(const ZFinder& zfind, AIDA::IHistogram1D* hist, BinnedHistogram<double>& binnedHist, double weight) {
        if (zfind.bosons().size() != 1) return;
        ParticleVector leptons = zfind.constituents();
        std::sort(leptons.begin(), leptons.end(), cmpParticleByPt);

        FourMomentum lminus = PID::threeCharge(leptons[0].pdgId()) < 0 ? leptons[0].momentum() : leptons[1].momentum();
        FourMomentum lplus = PID::threeCharge(leptons[0].pdgId()) < 0 ? leptons[1].momentum() : leptons[0].momentum();

        double phi_acop = M_PI - deltaPhi(lminus, lplus);
        double costhetastar = tanh((lminus.eta()-lplus.eta())/2.0);
        double sin2thetastar = 1.0 - sqr(costhetastar);
        if (sin2thetastar < 0.0) sin2thetastar = 0.0;
        double phistar = tan(phi_acop/2.0) * sqrt(sin2thetastar);
        hist->fill(phistar, weight);

        FourMomentum zmom = zfind.bosons()[0].momentum();
        binnedHist.fill(fabs(zmom.rapidity()), phistar, weight);
      }


      /// Normalise histograms etc., after the run
      void finalize() {
        normalize(_hist_zphistar_el_dressed);
        normalize(_hist_zphistar_el_bare);
        normalize(_hist_zphistar_mu_dressed);
        normalize(_hist_zphistar_mu_bare);

        foreach (AIDA::IHistogram1D* hist, _h_phistar_mu_dressed.getHistograms()) { normalize(hist); }
        foreach (AIDA::IHistogram1D* hist, _h_phistar_mu_bare.getHistograms()) { normalize(hist); }
        foreach (AIDA::IHistogram1D* hist, _h_phistar_el_bare.getHistograms()) { normalize(hist); }
        foreach (AIDA::IHistogram1D* hist, _h_phistar_el_dressed.getHistograms()) { normalize(hist); }
      }

      //@}


    private:

      BinnedHistogram<double> _h_phistar_mu_dressed;
      BinnedHistogram<double> _h_phistar_mu_bare;
      BinnedHistogram<double> _h_phistar_el_dressed;
      BinnedHistogram<double> _h_phistar_el_bare;

      AIDA::IHistogram1D *_hist_zphistar_el_dressed;
      AIDA::IHistogram1D *_hist_zphistar_el_bare;

      AIDA::IHistogram1D *_hist_zphistar_mu_dressed;
      AIDA::IHistogram1D *_hist_zphistar_mu_bare;

      AIDA::IHistogram1D *_hist_zphistar_el_bare_1;
      AIDA::IHistogram1D *_hist_zphistar_el_bare_2;
      AIDA::IHistogram1D *_hist_zphistar_el_bare_3;

      AIDA::IHistogram1D *_hist_zphistar_el_dressed_1;
      AIDA::IHistogram1D *_hist_zphistar_el_dressed_2;
      AIDA::IHistogram1D *_hist_zphistar_el_dressed_3;

      AIDA::IHistogram1D *_hist_zphistar_mu_bare_1;
      AIDA::IHistogram1D *_hist_zphistar_mu_bare_2;
      AIDA::IHistogram1D *_hist_zphistar_mu_bare_3;

      AIDA::IHistogram1D *_hist_zphistar_mu_dressed_1;
      AIDA::IHistogram1D *_hist_zphistar_mu_dressed_2;
      AIDA::IHistogram1D *_hist_zphistar_mu_dressed_3;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1204784);

}
