// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetYODA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class ATLAS_2011_S9128077 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    ATLAS_2011_S9128077()
      : Analysis("ATLAS_2011_S9128077")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;

      FastJets j4(fs, FastJets::ANTIKT, 0.4);
      j4.useInvisibles();
      addProjection(j4, "AntiKtJets04");

      FastJets j6(fs, FastJets::ANTIKT, 0.6);
      j6.useInvisibles();
      addProjection(j6, "AntiKtJets06");

      _h_jet_multi_inclusive = bookHisto1D(1, 1, 1);
      _h_jet_multi_ratio = bookScatter2D(2, 1, 1);
      _h_jet_pT.resize(4);
      _h_jet_pT[0] = bookHisto1D(3, 1, 1);
      _h_jet_pT[1] = bookHisto1D(4, 1, 1);
      _h_jet_pT[2] = bookHisto1D(5, 1, 1);
      _h_jet_pT[3] = bookHisto1D(6, 1, 1);
      _h_HT_2 = bookHisto1D(7, 1, 1);
      _h_HT_3 = bookHisto1D(8, 1, 1);
      _h_HT_4 = bookHisto1D(9, 1, 1);

      /// temporary histograms which will be divided in the end for the dsigma3/dsigma2 ratios
      _h_tmp_pTlead_R06_60_2  = bookHisto1D(10, 1, 1, "tmp1");
      _h_tmp_pTlead_R06_80_2  = bookHisto1D(11, 1, 1, "tmp2");
      _h_tmp_pTlead_R06_110_2 = bookHisto1D(12, 1, 1, "tmp3");
      _h_tmp_pTlead_R06_60_3  = bookHisto1D(10, 1, 1, "tmp4");
      _h_tmp_pTlead_R06_80_3  = bookHisto1D(11, 1, 1, "tmp5");
      _h_tmp_pTlead_R06_110_3 = bookHisto1D(12, 1, 1, "tmp6");

      _h_tmp_pTlead_R04_60_2  = bookHisto1D(13, 1, 1, "tmp7");
      _h_tmp_pTlead_R04_80_2  = bookHisto1D(14, 1, 1, "tmp8");
      _h_tmp_pTlead_R04_110_2 = bookHisto1D(15, 1, 1, "tmp9");
      _h_tmp_pTlead_R04_60_3  = bookHisto1D(13, 1, 1, "tmp10");
      _h_tmp_pTlead_R04_80_3  = bookHisto1D(14, 1, 1, "tmp11");
      _h_tmp_pTlead_R04_110_3 = bookHisto1D(15, 1, 1, "tmp12");

      _h_tmp_HT2_R06_2 = bookHisto1D(16, 1, 1, "tmp13");
      _h_tmp_HT2_R06_3 = bookHisto1D(16, 1, 1, "tmp14");
      _h_tmp_HT2_R04_2 = bookHisto1D(17, 1, 1, "tmp15");
      _h_tmp_HT2_R04_3 = bookHisto1D(17, 1, 1, "tmp16");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      vector<FourMomentum> jets04;
      foreach (const Jet& jet, applyProjection<FastJets>(event, "AntiKtJets04").jetsByPt(60.0*GeV)) {
        if (fabs(jet.momentum().eta())<2.8) {
          jets04.push_back(jet.momentum());
        }
      }

      if (jets04.size()>1 && jets04[0].pT()>80.0*GeV) {
        for (size_t i=2; i<=jets04.size(); ++i) {
          _h_jet_multi_inclusive->fill(i, weight);
        }

        double HT=0.0;
        for (size_t i=0; i<jets04.size(); ++i) {
          if (i<_h_jet_pT.size()) _h_jet_pT[i]->fill(jets04[i].pT(), weight);
          HT += jets04[i].pT();
        }

        if (jets04.size()>=2) _h_HT_2->fill(HT, weight);
        if (jets04.size()>=3) _h_HT_3->fill(HT, weight);
        if (jets04.size()>=4) _h_HT_4->fill(HT, weight);

        double pT1(jets04[0].pT()), pT2(jets04[1].pT());
        double HT2=pT1+pT2;
        if (jets04.size()>=2) {
          _h_tmp_HT2_R04_2->fill(HT2, weight);
          _h_tmp_pTlead_R04_60_2->fill(pT1, weight);
          if (pT2>80.0*GeV) _h_tmp_pTlead_R04_80_2->fill(pT1, weight);
          if (pT2>110.0*GeV) _h_tmp_pTlead_R04_110_2->fill(pT1, weight);
        }
        if (jets04.size()>=3) {
          double pT3(jets04[2].pT());
          _h_tmp_HT2_R04_3->fill(HT2, weight);
          _h_tmp_pTlead_R04_60_3->fill(pT1, weight);
          if (pT3>80.0*GeV) _h_tmp_pTlead_R04_80_3->fill(pT1, weight);
          if (pT3>110.0*GeV) _h_tmp_pTlead_R04_110_3->fill(pT1, weight);
        }
      }

      vector<FourMomentum> jets06;
      foreach (const Jet& jet, applyProjection<FastJets>(event, "AntiKtJets06").jetsByPt(60.0*GeV)) {
        if (fabs(jet.momentum().eta())<2.8) {
          jets06.push_back(jet.momentum());
        }
      }
      if (jets06.size()>1 && jets06[0].pT()>80.0*GeV) {
        double pT1(jets06[0].pT()), pT2(jets06[1].pT());
        double HT2=pT1+pT2;
        if (jets06.size()>=2) {
          _h_tmp_HT2_R06_2->fill(HT2, weight);
          _h_tmp_pTlead_R06_60_2->fill(pT1, weight);
          if (pT2>80.0*GeV) _h_tmp_pTlead_R06_80_2->fill(pT1, weight);
          if (pT2>110.0*GeV) _h_tmp_pTlead_R06_110_2->fill(pT1, weight);
        }
        if (jets06.size()>=3) {
          double pT3(jets06[2].pT());
          _h_tmp_HT2_R06_3->fill(HT2, weight);
          _h_tmp_pTlead_R06_60_3->fill(pT1, weight);
          if (pT3>80.0*GeV) _h_tmp_pTlead_R06_80_3->fill(pT1, weight);
          if (pT3>110.0*GeV) _h_tmp_pTlead_R06_110_3->fill(pT1, weight);
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Fill inclusive jet multi ratio
      int Nbins = _h_jet_multi_inclusive->numBins();
      std::vector<double> ratio(Nbins-1, 0.0);
      std::vector<double> err(Nbins-1, 0.0);
      for (int i = 0; i < Nbins-1; ++i) {
        if (_h_jet_multi_inclusive->bin(i).area() > 0.0 && _h_jet_multi_inclusive->bin(i+1).area() > 0.0) {
          ratio[i] = _h_jet_multi_inclusive->bin(i+1).area()/_h_jet_multi_inclusive->bin(i).area();
          double relerr_i = _h_jet_multi_inclusive->bin(i).areaErr()/_h_jet_multi_inclusive->bin(i).area();
          double relerr_j = _h_jet_multi_inclusive->bin(i+1).areaErr()/_h_jet_multi_inclusive->bin(i+1).area();
          err[i] = ratio[i] * (relerr_i + relerr_j);
        }
      }
      //\todo YODA
      //_h_jet_multi_ratio->setCoordinate(1, ratio, err);

      scale(_h_jet_multi_inclusive, crossSectionPerEvent());
      scale(_h_jet_pT[0], crossSectionPerEvent());
      scale(_h_jet_pT[1], crossSectionPerEvent());
      scale(_h_jet_pT[2], crossSectionPerEvent());
      scale(_h_jet_pT[3], crossSectionPerEvent());
      scale(_h_HT_2, crossSectionPerEvent());
      scale(_h_HT_3, crossSectionPerEvent());
      scale(_h_HT_4, crossSectionPerEvent());

      /// create ratio histograms
      //\todo YODA divide
      //histogramFactory().divide(histoDir() + "/d10-x01-y01", *_h_tmp_pTlead_R06_60_3, *_h_tmp_pTlead_R06_60_2);
      //histogramFactory().divide(histoDir() + "/d11-x01-y01", *_h_tmp_pTlead_R06_80_3, *_h_tmp_pTlead_R06_80_2);
      //histogramFactory().divide(histoDir() + "/d12-x01-y01", *_h_tmp_pTlead_R06_110_3, *_h_tmp_pTlead_R06_110_2);
      //histogramFactory().divide(histoDir() + "/d13-x01-y01", *_h_tmp_pTlead_R04_60_3, *_h_tmp_pTlead_R04_60_2);
      //histogramFactory().divide(histoDir() + "/d14-x01-y01", *_h_tmp_pTlead_R04_80_3, *_h_tmp_pTlead_R04_80_2);
      //histogramFactory().divide(histoDir() + "/d15-x01-y01", *_h_tmp_pTlead_R04_110_3, *_h_tmp_pTlead_R04_110_2);
      //histogramFactory().divide(histoDir() + "/d16-x01-y01", *_h_tmp_HT2_R06_3, *_h_tmp_HT2_R06_2);
      //histogramFactory().divide(histoDir() + "/d17-x01-y01", *_h_tmp_HT2_R04_3, *_h_tmp_HT2_R04_2);
      //histogramFactory().destroy(_h_tmp_pTlead_R06_60_2);
      //histogramFactory().destroy(_h_tmp_pTlead_R06_80_2);
      //histogramFactory().destroy(_h_tmp_pTlead_R06_110_2);
      //histogramFactory().destroy(_h_tmp_pTlead_R06_60_3);
      //histogramFactory().destroy(_h_tmp_pTlead_R06_80_3);
      //histogramFactory().destroy(_h_tmp_pTlead_R06_110_3);
      //histogramFactory().destroy(_h_tmp_pTlead_R04_60_2);
      //histogramFactory().destroy(_h_tmp_pTlead_R04_80_2);
      //histogramFactory().destroy(_h_tmp_pTlead_R04_110_2);
      //histogramFactory().destroy(_h_tmp_pTlead_R04_60_3);
      //histogramFactory().destroy(_h_tmp_pTlead_R04_80_3);
      //histogramFactory().destroy(_h_tmp_pTlead_R04_110_3);
      //histogramFactory().destroy(_h_tmp_HT2_R06_2);
      //histogramFactory().destroy(_h_tmp_HT2_R06_3);
      //histogramFactory().destroy(_h_tmp_HT2_R04_2);
      //histogramFactory().destroy(_h_tmp_HT2_R04_3);

    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_jet_multi_inclusive;
    Scatter2DPtr _h_jet_multi_ratio;
    vector<Histo1DPtr> _h_jet_pT;
    Histo1DPtr _h_HT_2;
    Histo1DPtr _h_HT_3;
    Histo1DPtr _h_HT_4;

    /// temporary histograms which will be divided in the end for the dsigma3/dsigma2 ratios
    Histo1DPtr _h_tmp_pTlead_R06_60_2;
    Histo1DPtr _h_tmp_pTlead_R06_80_2;
    Histo1DPtr _h_tmp_pTlead_R06_110_2;
    Histo1DPtr _h_tmp_pTlead_R06_60_3;
    Histo1DPtr _h_tmp_pTlead_R06_80_3;
    Histo1DPtr _h_tmp_pTlead_R06_110_3;

    Histo1DPtr _h_tmp_pTlead_R04_60_2;
    Histo1DPtr _h_tmp_pTlead_R04_80_2;
    Histo1DPtr _h_tmp_pTlead_R04_110_2;
    Histo1DPtr _h_tmp_pTlead_R04_60_3;
    Histo1DPtr _h_tmp_pTlead_R04_80_3;
    Histo1DPtr _h_tmp_pTlead_R04_110_3;

    Histo1DPtr _h_tmp_HT2_R06_2;
    Histo1DPtr _h_tmp_HT2_R06_3;
    Histo1DPtr _h_tmp_HT2_R04_2;
    Histo1DPtr _h_tmp_HT2_R04_3;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2011_S9128077);


}
