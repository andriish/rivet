// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Math/MathUtils.hh"
#include "LWH/Histogram1D.h"

namespace Rivet {


  namespace {

    inline double sumAB(vector<double> vecX, vector<double> vecY, vector<double> vecW) {
      assert(vecX.size() == vecY.size() && vecX.size() == vecW.size());
      double sum(0);
      for (size_t i = 0; i < vecX.size(); i++) sum += vecW[i] * vecX[i] * vecY[i];
      return sum;
    }

    inline double sumA(vector<double> vecX, vector<double> vecW) {
      assert(vecX.size() == vecW.size());
      double sum(0);
      for (size_t i = 0; i < vecX.size(); i++) sum += vecX[i]*vecW[i];
      return sum;
    }

    inline double sumW(vector<double> vecW) {
      double sum(0);
      for (size_t i = 0; i < vecW.size(); i++) sum += vecW[i];
      return sum;
    }

    inline double mean(vector<double>  vecX, vector<double> vecW) {
      return sumA(vecX, vecW) / sumW(vecW);
    }

    inline double standard_deviation(vector<double> vecX, vector<double> vecW) {
      const double x_bar = mean(vecX, vecW);
      double sum(0);
      for (size_t i = 0; i < vecX.size(); i++) {
        sum += vecW[i] * sqr(vecX[i] - x_bar);
      }
      return sqrt( sum / sumW(vecW) );
    }

    inline double a0_regression(vector<double>  vecX, vector<double> vecY, vector<double> vecW) {
      const double numerator   = sumA(vecY, vecW) * sumAB(vecX, vecX, vecW) - sumA(vecX, vecW) * sumAB(vecX, vecY, vecW);
      const double denominator = sumW(vecW) * sumAB(vecX, vecX, vecW) - sumA(vecX, vecW) * sumA(vecX, vecW);
      return numerator / denominator;
    }

    inline double a1_regression(vector<double>  vecX, vector<double> vecY, vector<double> vecW) {
      const double numerator   = sumW(vecW) * sumAB(vecX,vecY,vecW) - sumA(vecX, vecW) * sumA(vecY, vecW);
      const double denominator = sumW(vecW) * sumAB(vecX,vecX,vecW) - sumA(vecX, vecW) * sumA(vecX, vecW);
      return numerator/ denominator;
    }

    inline double a1_regression2(vector<double>  vecX, vector<double> vecY, vector<double> vecW) {
      const double x_bar = mean(vecX, vecW);
      const double y_bar = mean(vecY, vecW);
      double sumXY(0);
      for (size_t i = 0; i < vecX.size(); i++) {
        sumXY += vecW[i] * (vecY[i]-y_bar) * (vecX[i]-x_bar);
      }
      return sumXY / ( standard_deviation(vecX, vecW) * standard_deviation(vecY, vecW) * sumW(vecW) );
    }

    inline double quadra_sum_residual(vector<double> vecX, vector<double> vecY, vector<double> vecW) {
      const double a0 = a0_regression(vecX, vecY, vecW);
      const double a1 = a1_regression(vecX, vecY, vecW);
      double sum(0);
      for (size_t i = 0; i < vecX.size(); i++) {
        const double y_est = a0 + a1*vecX[i];
        sum += vecW[i] * sqr(vecY[i] - y_est);
      }
      return sum;
    }

    inline double error_on_slope(vector<double> vecX, vector<double> vecY, vector<double> vecW) {
      const double quadra_sum_res = quadra_sum_residual(vecX, vecY, vecW);
      const double sqrt_quadra_sum_x   = standard_deviation(vecX, vecW) * sqrt(sumW(vecW));
      return sqrt(quadra_sum_res/(sumW(vecW)-2)) / sqrt_quadra_sum_x;
    }

  }



  /// Generic analysis looking at various distributions of final state particles
  class ATLAS_2012_I1093734 : public Analysis {
  public:

    /// Constructor
    ATLAS_2012_I1093734()
      : Analysis("ATLAS_2012_I1093734")
    {
      // Stat convergence happens around 20k events, so it doesn't make sense to run this
      // analysis with much less than that. Given that, lets avoid some unnecessary vector
      // resizing by allocating sensible amounts in the first place.
      for (int ipt = 0; ipt < NptBins; ++ipt) {
        for (int k = 0; k < NetaBins; ++k) {
          _vecsNchF [ipt][k].reserve(10000);
          _vecsNchB [ipt][k].reserve(10000);
          _vecWeight[ipt][k].reserve(10000);
          if (ipt == 0) {
            _vecsSumptF[k].reserve(10000);
            _vecsSumptB[k].reserve(10000);
          }
        }
      }
    }


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // FB correlations part

      // Projections
      for (int ipt = 0; ipt < NptBins; ++ipt) {
        const double ptmin = ptminValue[ipt]*MeV;
        for (int ieta = 0; ieta < NetaBins; ++ieta) {
          addProjection(ChargedFinalState(-etaValue[ieta]    , -etaValue[ieta]+0.5, ptmin), "Tracks"+etabinnames[ieta]+"B"+ptbinnames[ipt]);
          addProjection(ChargedFinalState( etaValue[ieta]-0.5,  etaValue[ieta]    , ptmin), "Tracks"+etabinnames[ieta]+"F"+ptbinnames[ipt]);
        }
        addProjection(ChargedFinalState(-2.5, 2.5, ptmin), "CFS" + ptbinnames[ipt]);
      }
      // Histos
      if (fuzzyEquals(sqrtS(), 7000*GeV, 1e-3)) {
        for (int ipt  = 0; ipt  < NptBins ; ++ipt )  _histNchCorr_vsEta[ipt ]= bookDataPointSet( 1+ipt ,2,1);
        for (int ieta = 0; ieta < NetaBins; ++ieta)  _histNchCorr_vsPt [ieta]= bookDataPointSet( 8+ieta,2,1);
        _histPtsumCorr = bookDataPointSet(13,2,1);
      } else if (fuzzyEquals(sqrtS(), 900*GeV, 1e-3)) {
        _histNchCorr_vsEta[0] = bookDataPointSet(14, 2, 1);
        _histPtsumCorr        = bookDataPointSet(15, 2, 1);
      }


      // Azimuthal correlations part
      // Projections
      const double ptmin = 500*MeV;
      addProjection(ChargedFinalState(-2.5, 2.5, ptmin), "ChargedTracks25");
      addProjection(ChargedFinalState(-2.0, 2.0, ptmin), "ChargedTracks20");
      addProjection(ChargedFinalState(-1.0, 1.0, ptmin), "ChargedTracks10");
      // Histos
      h_dphi[0] = bookHistogram1D("dphi0", 50, 0.0, PI);
      h_same[0] = bookHistogram1D("same0", 50, 0.0, PI);
      h_oppo[0] = bookHistogram1D("oppo0", 50, 0.0, PI);
      h_dphi[1] = bookHistogram1D("dphi1", 50, 0.0, PI);
      h_same[1] = bookHistogram1D("same1", 50, 0.0, PI);
      h_oppo[1] = bookHistogram1D("oppo1", 50, 0.0, PI);
      h_dphi[2] = bookHistogram1D("dphi2", 50, 0.0, PI);
      h_same[2] = bookHistogram1D("same2", 50, 0.0, PI);
      h_oppo[2] = bookHistogram1D("oppo2", 50, 0.0, PI);
      if (fuzzyEquals(sqrtS(), 7000*GeV, 1e-3)) {
        h_dphiMin[0] = bookDataPointSet(2,1,1);
        h_dphiMin[1] = bookDataPointSet(4,1,1);
        h_dphiMin[2] = bookDataPointSet(6,1,1);
      } else if (fuzzyEquals(sqrtS(), 900*GeV, 1e-3)) {
        h_dphiMin[0] = bookDataPointSet(1,1,1);
        h_dphiMin[1] = bookDataPointSet(3,1,1);
        h_dphiMin[2] = bookDataPointSet(5,1,1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      for (int ipt = 0; ipt < NptBins; ++ipt) {
        const FinalState& charged = applyProjection<FinalState>(event, "CFS" + ptbinnames[ipt]);
        if (charged.particles().size() >= 2) {
          for (int ieta = 0; ieta < NetaBins; ++ieta) {
            const string fname = "Tracks" + etabinnames[ieta] + "F" + ptbinnames[ipt];
            const string bname = "Tracks" + etabinnames[ieta] + "B" + ptbinnames[ipt];
            const ParticleVector particlesF = applyProjection<FinalState>(event, fname).particles();
            const ParticleVector particlesB = applyProjection<FinalState>(event, bname).particles();
            _vecsNchF[ipt][ieta].push_back((double) particlesF.size());
            _vecsNchB[ipt][ieta].push_back((double) particlesB.size());
            _vecWeight[ipt][ieta].push_back(weight);

            // Sum pT only for 100 MeV particles
            if (ipt == 0) {
              double sumptF = 0;
              double sumptB = 0;
              foreach (const Particle& p, particlesF) sumptF += p.momentum().pT();
              foreach (const Particle& p, particlesB) sumptB += p.momentum().pT();
              _vecsSumptF[ieta].push_back(sumptF);
              _vecsSumptB[ieta].push_back(sumptB);
            }

          }
        }
      }


      string etabin[3] = { "10", "20", "25" };
      for (int iColl = 0; iColl < 3; iColl++) {
        const string fname = "ChargedTracks" + etabin[iColl];
        const ParticleVector partTrks = applyProjection<FinalState>(event, fname).particlesByPt();

        // First loop over particles to find the leading track
        const Particle& plead = partTrks[0];

        // Second loop to fill the histograms
        foreach (const Particle& p, partTrks) {
          if (&plead == &p) continue; //< Don't compare the lead particle to itself
          const double dphi = deltaPhi(p.momentum(), plead.momentum());
          h_dphi[iColl]->fill(dphi, weight);
          const bool sameside = (plead.momentum().eta() * p.momentum().eta() > 0);
          (sameside ? h_same : h_oppo)[iColl]->fill(dphi, weight);
        }
      }
    }


    /// Finalize
    void finalize() {

      // FB part
      // @todo For 2D plots we will need _vecsNchF[i], _vecsNchB[j]
      for (int ipt = 0; ipt < NptBins; ++ipt) {
        vector<double> nchval, ncherr;
        for (int k = 0; k < NetaBins; ++k) {
          nchval.push_back( a1_regression2(_vecsNchF[ipt][k], _vecsNchB[ipt][k], _vecWeight[ipt][k]));
          ncherr.push_back( error_on_slope(_vecsNchF[ipt][k], _vecsNchB[ipt][k], _vecWeight[ipt][k]));
        }
        _histNchCorr_vsEta[ipt]->setCoordinate(1, nchval, ncherr);
        // There is just one plot at 900 GeV so exit the loop here
        if (fuzzyEquals(sqrtS(),  900*GeV, 1e-3) && ipt > 0) break;
      }

      if (!fuzzyEquals(sqrtS(),  900*GeV, 1e-3)) { //< No plots at 900 GeV
        for (int ieta = 0; ieta < NetaBins; ++ieta) {
          vector<double> nchval, ncherr;
          for (int ipt = 0; ipt < NptBins; ++ipt) {
            nchval.push_back( a1_regression2(_vecsNchF[ipt][ieta], _vecsNchB[ipt][ieta], _vecWeight[ipt][ieta]));
            ncherr.push_back( error_on_slope(_vecsNchF[ipt][ieta], _vecsNchB[ipt][ieta], _vecWeight[ipt][ieta]));
          }
          _histNchCorr_vsPt [ieta]->setCoordinate(1, nchval, ncherr);
        }
      }

      // Sum pt only for 100 MeV particles
      vector<double> ptsval, ptserr;
      for (int k = 0; k < NetaBins; ++k) {
        ptsval.push_back(a1_regression2(_vecsSumptF[k], _vecsSumptB[k], _vecWeight[0][k]));
        ptserr.push_back(error_on_slope(_vecsSumptF[k], _vecsSumptB[k], _vecWeight[0][k]));
      }
      _histPtsumCorr->setCoordinate(1, ptsval, ptserr);


      // Azimuthal part
      AIDA::IHistogramFactory& hf = histogramFactory();
      AIDA::IHistogram1D* h_sameMinusOppo[3];
      for (int iColl = 0; iColl < 3; iColl++) {
        if (fuzzyEquals(sqrtS(), 7000*GeV, 1e-3)) {
          h_sameMinusOppo[iColl] = hf.subtract( histoPath(8+2*iColl,1,1), *h_same[iColl], *h_oppo[iColl]);
          scale(h_sameMinusOppo[iColl], 1./h_sameMinusOppo[iColl]->sumBinHeights() * h_sameMinusOppo[iColl]->axis().binWidth(1));
        } else if (fuzzyEquals(sqrtS(), 900*GeV, 1e-3)) {
          h_sameMinusOppo[iColl] = hf.subtract( histoPath(7+2*iColl,1,1), *h_same[iColl], *h_oppo[iColl]);
          scale(h_sameMinusOppo[iColl], 1./h_sameMinusOppo[iColl]->sumBinHeights() * h_sameMinusOppo[iColl]->axis().binWidth(1));
        }
        int Nbins      = h_dphi[iColl]->axis().bins();
        double histMin = h_dphi[iColl]->minBinHeight();

        std::vector<double> diff(Nbins, 0.0);
        std::vector<double> err(Nbins, 0.0);

        // Subtract minimal value
        double sumDiff = 0.;
        for (int n = 0; n < Nbins; ++n) {
          diff[n] = h_dphi[iColl]->binHeight(n) - histMin;
          err[n]  = h_dphi[iColl]->binError(n);
          sumDiff += diff[n];
        }
        // Normalize
        for (int n = 0; n < Nbins; ++n) {
          diff[n] = diff[n]/sumDiff;
          err[n]  = err[n]/sumDiff;
        }
        h_dphiMin[iColl]->setCoordinate(1, diff, err);
      }

      // Remove temporary histograms
      for (size_t i = 0; i < 3; i++) {
        hf.destroy(h_dphi[i]);
        hf.destroy(h_same[i]);
        hf.destroy(h_oppo[i]);
      }

    }

    //@}


  private:

    static const int NptBins  = 7;
    static const int NetaBins = 5;
    static const double ptminValue[NptBins];
    static const string ptbinnames[NptBins];
    static const double etaValue   [NetaBins];
    static const string etabinnames[NetaBins];

    vector<double> _vecWeight[NptBins][NetaBins];
    vector<double> _vecsNchF [NptBins][NetaBins];
    vector<double> _vecsNchB [NptBins][NetaBins];
    vector<double> _vecsSumptF[NetaBins];
    vector<double> _vecsSumptB[NetaBins];

    /// @name Histograms
    //@{
    AIDA::IDataPointSet *_histNchCorr_vsEta[NptBins ];
    AIDA::IDataPointSet *_histNchCorr_vsPt [NetaBins];
    AIDA::IDataPointSet *_histPtsumCorr;
    AIDA::IHistogram1D* h_dphi[3];
    AIDA::IHistogram1D* h_same[3];
    AIDA::IHistogram1D* h_oppo[3];
    AIDA::IDataPointSet* h_dphiMin[3];
    //@}

  };


  /// @todo Initialize these inline at declaration with C++11
  const double ATLAS_2012_I1093734::ptminValue[] = {100, 200, 300, 500, 1000, 1500, 2000 };
  const string ATLAS_2012_I1093734::ptbinnames[] = { "100", "200", "300", "500", "1000", "1500", "2000" };
  const double ATLAS_2012_I1093734::etaValue   [] = {0.5, 1.0, 1.5, 2.0, 2.5};
  const string ATLAS_2012_I1093734::etabinnames[] = { "05", "10", "15", "20", "25" };



  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1093734);

}
