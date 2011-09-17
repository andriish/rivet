// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"


namespace Rivet {


  class D0_2010_S8821313 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    D0_2010_S8821313()
      : Analysis("D0_2010_S8821313")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// Initialise and register projections
      vector<pair<double, double> > etaRanges_ee;
      etaRanges_ee.push_back(make_pair(-3.0, -1.5));
      etaRanges_ee.push_back(make_pair(-1.1, 1.1));
      etaRanges_ee.push_back(make_pair(1.5, 3.0));
      ZFinder zfinder_ee(etaRanges_ee, 20.0*GeV, ELECTRON, 70.0*GeV, 110.0*GeV, 0.2, true, true);
      addProjection(zfinder_ee, "zfinder_ee");

      ZFinder zfinder_mm(-2.0, 2.0, 15.0*GeV, MUON, 70.0*GeV, 110.0*GeV, 0.0, false, false);
      addProjection(zfinder_mm, "zfinder_mm");

      /// Book histograms here
      _h_phistar_ee.addHistogram(0.0, 1.0, bookHistogram1D(1, 1, 1));
      _h_phistar_ee.addHistogram(1.0, 2.0, bookHistogram1D(1, 1, 2));
      _h_phistar_ee.addHistogram(2.0, 10.0, bookHistogram1D(1, 1, 3));

      _h_phistar_mm.addHistogram(0.0, 1.0, bookHistogram1D(2, 1, 1));
      _h_phistar_mm.addHistogram(1.0, 2.0, bookHistogram1D(2, 1, 2));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      const ZFinder& zfinder_ee = applyProjection<ZFinder>(event, "zfinder_ee");
      if (zfinder_ee.bosons().size()==1) {
        ParticleVector ee=zfinder_ee.constituents();
        std::sort(ee.begin(), ee.end(), cmpParticleByPt);
        FourMomentum eminus=PID::threeCharge(ee[0].pdgId())<0.0?ee[0].momentum():ee[1].momentum();
        FourMomentum eplus=PID::threeCharge(ee[0].pdgId())<0.0?ee[1].momentum():ee[0].momentum();
        double phi_acop=M_PI-mapAngle0ToPi(eminus.phi()-eplus.phi());
        double costhetastar=tanh((eminus.eta()-eplus.eta())/2.0);
        double sin2thetastar=1.0-sqr(costhetastar);
        if (sin2thetastar<0.0) sin2thetastar=0.0;
        double phistar=tan(phi_acop/2.0)*sqrt(sin2thetastar);

        FourMomentum Zmom=zfinder_ee.bosons()[0].momentum();
        _h_phistar_ee.fill(Zmom.rapidity(), phistar, weight);
      }

      const ZFinder& zfinder_mm = applyProjection<ZFinder>(event, "zfinder_mm");
      if (zfinder_mm.bosons().size()==1) {
        ParticleVector mm=zfinder_mm.constituents();
        std::sort(mm.begin(), mm.end(), cmpParticleByPt);
        FourMomentum mminus=PID::threeCharge(mm[0].pdgId())<0.0?mm[0].momentum():mm[1].momentum();
        FourMomentum mplus=PID::threeCharge(mm[0].pdgId())<0.0?mm[1].momentum():mm[0].momentum();
        double phi_acop=M_PI-mapAngle0ToPi(mminus.phi()-mplus.phi());
        double costhetastar=tanh((mminus.eta()-mplus.eta())/2.0);
        double sin2thetastar=1.0-sqr(costhetastar);
        if (sin2thetastar<0.0) sin2thetastar=0.0;
        double phistar=tan(phi_acop/2.0)*sqrt(sin2thetastar);

        FourMomentum Zmom=zfinder_mm.bosons()[0].momentum();
        _h_phistar_mm.fill(Zmom.rapidity(), phistar, weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      foreach (AIDA::IHistogram1D* hist, _h_phistar_ee.getHistograms()) {
        normalize(hist, 1.0);
      }
      foreach (AIDA::IHistogram1D* hist, _h_phistar_mm.getHistograms()) {
        normalize(hist, 1.0);
      }
    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here


  private:

    /// @name Histograms
    //@{

    BinnedHistogram<double> _h_phistar_ee;
    BinnedHistogram<double> _h_phistar_mm;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(D0_2010_S8821313);

}
