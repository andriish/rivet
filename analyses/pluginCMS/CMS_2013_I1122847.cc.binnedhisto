// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  class CMS_2013_I1122847 : public Analysis {
  public:

    /// Constructor
    CMS_2013_I1122847()
      : Analysis("CMS_2013_I1122847")  {}


    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs;

      Cut cuts_mu = Cuts::abseta < 2.4 && Cuts::pT > 20*GeV;
      ZFinder zfinder_mu(fs, cuts_mu, PID::MUON, 40.0*GeV, MAXDOUBLE,
                         0.0, ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);
      declare(zfinder_mu, "zfinder_mu");

      Cut cuts_el = (Cuts::pT >= 20*GeV && Cuts::abseta < 2.4 && !Cuts::absetaIn(1.447, 1.57));
      ZFinder zfinder_el(fs, cuts_el, PID::ELECTRON, 40.0*GeV, MAXDOUBLE,
                         0.0, ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);
      declare(zfinder_el, "zfinder_el");


      /// Histograms
      // dimuon
      book(_hist_mm_100_num, "TMP/mm_100_num", refData(1, 1, 1));
      book(_hist_mm_125_num, "TMP/mm_125_num", refData(1, 1, 2));
      book(_hist_mm_150_num, "TMP/mm_150_num", refData(1, 1, 3));
      book(_hist_mm_240_num, "TMP/mm_240_num", refData(1, 1, 4));

      book(_hist_mm_100_den, "TMP/mm_100_den", refData(1, 1, 1));
      book(_hist_mm_125_den, "TMP/mm_125_den", refData(1, 1, 2));
      book(_hist_mm_150_den, "TMP/mm_150_den", refData(1, 1, 3));
      book(_hist_mm_240_den, "TMP/mm_240_den", refData(1, 1, 4));

      // Dielectron
      book(_hist_ee_100_num, "TMP/ee_100_num", refData(2, 1, 1));
      book(_hist_ee_125_num, "TMP/ee_125_num", refData(2, 1, 2));
      book(_hist_ee_150_num, "TMP/ee_150_num", refData(2, 1, 3));
      book(_hist_ee_240_num, "TMP/ee_240_num", refData(2, 1, 4));

      book(_hist_ee_100_den, "TMP/ee_100_den", refData(2, 1, 1));
      book(_hist_ee_125_den, "TMP/ee_125_den", refData(2, 1, 2));
      book(_hist_ee_150_den, "TMP/ee_150_den", refData(2, 1, 3));
      book(_hist_ee_240_den, "TMP/ee_240_den", refData(2, 1, 4));

      // Dilepton
      book(_hist_ll_100_num, "TMP/ll_100_num", refData(3, 1, 1));
      book(_hist_ll_125_num, "TMP/ll_125_num", refData(3, 1, 2));
      book(_hist_ll_150_num, "TMP/ll_150_num", refData(3, 1, 3));
      book(_hist_ll_240_num, "TMP/ll_240_num", refData(3, 1, 4));

      book(_hist_ll_100_den, "TMP/ll_100_den", refData(3, 1, 1));
      book(_hist_ll_125_den, "TMP/ll_125_den", refData(3, 1, 2));
      book(_hist_ll_150_den, "TMP/ll_150_den", refData(3, 1, 3));
      book(_hist_ll_240_den, "TMP/ll_240_den", refData(3, 1, 4));
    }


    double cosThetaCS(const Particle& l1, const Particle& l2) {
      const FourMomentum mom1 = l1.mom();
      const FourMomentum mom2 = l2.mom();
      const FourMomentum mom12 = mom1 + mom2;
      const double Q = mom12.mass();
      const double QT = mom12.pT();
      const double QZ = mom12.pz();

      /// @todo Why include factors of sqrt2 which then get immediately multiplied then divided out?
      const double sqrt2 = sqrt(2.0);
      /// @todo Can be done more nicely via PID-ordered references to mom1, mom2
      const double P1p = ((l1.pid() > 0) ? (mom1.E() + mom1.pz()) : (mom2.E() + mom2.pz())) / sqrt2;
      const double P1m = ((l1.pid() > 0) ? (mom1.E() - mom1.pz()) : (mom2.E() - mom2.pz())) / sqrt2;
      const double P2p = ((l1.pid() > 0) ? (mom2.E() + mom2.pz()) : (mom1.E() + mom1.pz())) / sqrt2;
      const double P2m = ((l1.pid() > 0) ? (mom2.E() - mom2.pz()) : (mom1.E() - mom1.pz())) / sqrt2;

      const double cosThetaCS = sign(QZ) * (2 / (Q * add_quad(Q, QT))) * (P1p*P2m - P1m*P2p);
      return cosThetaCS;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      const ZFinder& zfinder_el = apply<ZFinder>(event, "zfinder_el");
      if (zfinder_el.bosons().size() > 0) {
        const Particle& z  = zfinder_el.bosons()[0];
        const Particle& l1 = zfinder_el.constituents()[0];
        const Particle& l2 = zfinder_el.constituents()[1];

        // Prepare variables for filling
        const double rap = z.absrap();
        const double costhetacs = cosThetaCS(l1, l2);
        const double sgn = sign(costhetacs);

        // Fill the histograms
        if (rap < 1.0) {
          _hist_ee_100_num->fill(z.mass(), weight * sgn);
          _hist_ll_100_num->fill(z.mass(), weight * sgn);
          _hist_ee_100_den->fill(z.mass(), weight);
          _hist_ll_100_den->fill(z.mass(), weight);
        } else if (rap < 1.25) {
          _hist_ee_125_num->fill(z.mass(), weight * sgn);
          _hist_ll_125_num->fill(z.mass(), weight * sgn);
          _hist_ee_125_den->fill(z.mass(), weight);
          _hist_ll_125_den->fill(z.mass(), weight);
        } else if (rap < 1.50) {
          _hist_ee_150_num->fill(z.mass(), weight * sgn);
          _hist_ll_150_num->fill(z.mass(), weight * sgn);
          _hist_ee_150_den->fill(z.mass(), weight);
          _hist_ll_150_den->fill(z.mass(), weight);
        } else if (rap < 2.40) {
          _hist_ee_240_num->fill(z.mass(), weight * sgn);
          _hist_ll_240_num->fill(z.mass(), weight * sgn);
          _hist_ee_240_den->fill(z.mass(), weight);
          _hist_ll_240_den->fill(z.mass(), weight);
        }
      }

      const ZFinder& zfinder_mu = apply<ZFinder>(event, "zfinder_mu");
      if (zfinder_mu.bosons().size() > 0) {
        const Particle& z  = zfinder_mu.bosons()[0];
        const Particle& l1 = zfinder_mu.constituents()[0];
        const Particle& l2 = zfinder_mu.constituents()[1];

        // Prepare variables for filling
        const double rap = z.absrap();
        const double costhetacs = cosThetaCS(l1, l2);
        const double sgn = sign(costhetacs);

        // Fill the histograms
        if (rap < 1.0) {
          _hist_mm_100_num->fill(z.mass(), weight * sgn);
          _hist_ll_100_num->fill(z.mass(), weight * sgn);
          _hist_mm_100_den->fill(z.mass(), weight);
          _hist_ll_100_den->fill(z.mass(), weight);
        } else if (rap < 1.25) {
          _hist_mm_125_num->fill(z.mass(), weight * sgn);
          _hist_ll_125_num->fill(z.mass(), weight * sgn);
          _hist_mm_125_den->fill(z.mass(), weight);
          _hist_ll_125_den->fill(z.mass(), weight);
        } else if (rap < 1.50) {
          _hist_mm_150_num->fill(z.mass(), weight * sgn);
          _hist_ll_150_num->fill(z.mass(), weight * sgn);
          _hist_mm_150_den->fill(z.mass(), weight);
          _hist_ll_150_den->fill(z.mass(), weight);
        } else if (rap < 2.40) {
          _hist_mm_240_num->fill(z.mass(), weight * sgn);
          _hist_ll_240_num->fill(z.mass(), weight * sgn);
          _hist_mm_240_den->fill(z.mass(), weight);
          _hist_ll_240_den->fill(z.mass(), weight);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      {Scatter2DPtr s2d; divide(_hist_mm_100_num, _hist_mm_100_den, book(s2d, 1, 1, 1));}
      {Scatter2DPtr s2d; divide(_hist_mm_125_num, _hist_mm_125_den, book(s2d, 1, 1, 2));}
      {Scatter2DPtr s2d; divide(_hist_mm_150_num, _hist_mm_150_den, book(s2d, 1, 1, 3));}
      {Scatter2DPtr s2d; divide(_hist_mm_240_num, _hist_mm_240_den, book(s2d, 1, 1, 4));}
      {Scatter2DPtr s2d; divide(_hist_ee_100_num, _hist_ee_100_den, book(s2d, 2, 1, 1));}
      {Scatter2DPtr s2d; divide(_hist_ee_125_num, _hist_ee_125_den, book(s2d, 2, 1, 2));}
      {Scatter2DPtr s2d; divide(_hist_ee_150_num, _hist_ee_150_den, book(s2d, 2, 1, 3));}
      {Scatter2DPtr s2d; divide(_hist_ee_240_num, _hist_ee_240_den, book(s2d, 2, 1, 4));}
      {Scatter2DPtr s2d; divide(_hist_ll_100_num, _hist_ll_100_den, book(s2d, 3, 1, 1));}
      {Scatter2DPtr s2d; divide(_hist_ll_125_num, _hist_ll_125_den, book(s2d, 3, 1, 2));}
      {Scatter2DPtr s2d; divide(_hist_ll_150_num, _hist_ll_150_den, book(s2d, 3, 1, 3));}
      {Scatter2DPtr s2d; divide(_hist_ll_240_num, _hist_ll_240_den, book(s2d, 3, 1, 4));}
    }


  private:

    /// Histograms
    Histo1DPtr _hist_ee_100_num, _hist_ee_125_num, _hist_ee_150_num, _hist_ee_240_num;
    Histo1DPtr _hist_ee_100_den, _hist_ee_125_den, _hist_ee_150_den, _hist_ee_240_den;
    Histo1DPtr _hist_mm_100_num, _hist_mm_125_num, _hist_mm_150_num, _hist_mm_240_num;
    Histo1DPtr _hist_mm_100_den, _hist_mm_125_den, _hist_mm_150_den, _hist_mm_240_den;
    Histo1DPtr _hist_ll_100_num, _hist_ll_125_num, _hist_ll_150_num, _hist_ll_240_num;
    Histo1DPtr _hist_ll_100_den, _hist_ll_125_den, _hist_ll_150_den, _hist_ll_240_den;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1122847);

}
