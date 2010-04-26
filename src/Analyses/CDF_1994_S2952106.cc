// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief CDF Run I color coherence analysis
  /// @author Lars Sonnenschein
  /// @todo BROKEN!
  class CDF_1994_S2952106 : public Analysis {
  public:

    /// Constructor
    CDF_1994_S2952106() : Analysis("CDF_1994_S2952106")
    {
      setBeams(PROTON, ANTIPROTON);
      // setNeedsCrossSection(true);
    }


    /// @name Analysis methods
    //@{

    void init() {
      const FinalState fs(-4.2, 4.2);
      addProjection(fs, "FS");
      addProjection(FastJets(fs, FastJets::CDFJETCLU, 0.7), "ConeJets");
      addProjection(TotalVisibleMomentum(fs), "CalMET");

      // Veto (anti)neutrinos, and muons with pT above 1.0 GeV
      VetoedFinalState vfs(fs);
      vfs.vetoNeutrinos();
      vfs.addVetoPairDetail(MUON, 1.0*GeV, MAXDOUBLE);
      addProjection(vfs, "VFS");

      /// @todo Use histogram auto-booking

      //const string hname = "HvsDphi";
      //const string htitle = "H vs Delta phi";
      //_histHvsDphi = bookHistogram2D(hname, htitle, 40, -4., 4., 32, 0., 3.2);

      //const string hname2 = "RvsAlpha";
      //const string htitle2 = "R vs alpha";
      //_histRvsAlpha = bookHistogram2D(hname2, htitle2, 50, 0., 5., 32, -1.6, 1.6);

      //_histJet1Et  = bookHistogram1D("Jet1Et", 40, 0., 500.);
      _histJet1Et  = bookHistogram1D(1,1,1);
      //_histJet2Et  = bookHistogram1D("Jet2Et", 40, 0., 500.);
      _histJet2Et  = bookHistogram1D(2,1,1);
      //_histJet3eta = bookHistogram1D("Jet3eta", 42, -4., 4.);
      _histJet3eta = bookHistogram1D(3,1,1);
      //_histR23     = bookHistogram1D("R23", 50, 0., 5.);
      _histR23     = bookHistogram1D(4,1,1);

      /// @todo Need better title
      //_histAlpha = bookHistogram1D("alpha", 42, -PI/2., PI/2.);
      _histAlpha = bookHistogram1D(5,1,1);

      //const string hname8 = "alphaMCvsDat";
      //const string htitle8 = "alpha MC vs. Data ";
      //_histAlphaMCvsDat = bookHistogram1D(hname8, htitle8, 42, -PI/2., PI/2.);

      /// @todo Need better title
      //_histAlphaCDF = bookHistogram1D("alphaCDF", 42, -PI/2., PI/2.);
      _histAlphaIdeal = bookHistogram1D(6,1,1);

      /// @todo Need better title
      //_histAlphaIdeal = bookHistogram1D("alphaIdeal", 42, -PI/2., PI/2.);
      _histAlphaIdeal = bookHistogram1D(6,1,2);


      /// @todo Need better title
      //_histR23CDF = bookHistogram1D("R23CDF", 50, 0., 5.);
      _histR23Ideal = bookHistogram1D(7,1,1);

      /// @todo Need better title
      //_histR23Ideal = bookHistogram1D("R23Ideal", 50, 0., 5.);
      _histR23Ideal = bookHistogram1D(7,1,2);
      

      /// @todo Need better title
      //_histJet3etaCDF = bookHistogram1D("Jet3etaCDF", 42, -4., 4.);
      _histJet3etaCDF = bookHistogram1D(8,1,1);

      /// @todo Need better title
      //_histJet3etaIdeal = bookHistogram1D("Jet3etaIdeal", 42, -4., 4.);
      _histJet3etaIdeal = bookHistogram1D(8,1,2);

    }



    // Do the analysis
    void analyze(const Event & event) {
      const Jets jets = applyProjection<FastJets>(event, "ConeJets").jetsByPt();
      getLog() << Log::DEBUG << "Jet multiplicity before any cuts = " << jets.size() << endl;

      // Check there isn't too much missing Et
      const TotalVisibleMomentum& caloMissEt = applyProjection<TotalVisibleMomentum>(event, "CalMET");
      getLog() << Log::DEBUG << "Missing pT = " << caloMissEt.momentum().pT()/GeV << " GeV" << endl;
      /// @todo should this really be scalar ET here, and not caloMissEt.momentum().Et()?
      //if ((caloMissEt.momentum().pT()/GeV)/sqrt(caloMissEt.scalarET()/GeV) > 6.0) vetoEvent;
      if ((caloMissEt.momentum().pT()/GeV)/sqrt(caloMissEt.momentum().Et()/GeV) > 6.0) vetoEvent;

      // Check jet requirements
      if (jets.size() < 3) vetoEvent;
      if (jets[0].momentum().pT() < 100*GeV) vetoEvent;

      // More jet 1,2,3 checks
      FourMomentum pj1(jets[0].momentum()), pj2(jets[1].momentum()), pj3(jets[2].momentum());
      if (fabs(pj1.eta()) > 0.7 || fabs(pj2.eta()) > 0.7) vetoEvent;
      getLog() << Log::DEBUG << "Jet 1 & 2 eta, pT requirements fulfilled" << endl;

      // back-to-bak within 20 degrees in phi
      if (deltaPhi(pj1.phi(), pj2.phi()) < PI-PI/18.) vetoEvent;
      getLog() << Log::DEBUG << "Jet 1 & 2 phi requirement fulfilled" << endl;

      const double weight = event.weight();
      _histJet1Et->fill(pj1.pT(), weight);
      _histJet2Et->fill(pj2.pT(), weight);
      _histR23->fill(deltaR(pj2, pj3), weight);
      _histJet3eta->fill(pj3.eta(), weight);

      // Next cut only required for alpha studies
      if (pj3.pT()/GeV < 10.0) vetoEvent;
      getLog() << Log::DEBUG << "3rd jet passes alpha histo pT cut" << endl;

      // Calc and plot alpha
      const double dPhi = deltaPhi(pj3.phi(), pj2.phi());
      const double dH = sign(pj2.eta()) * (pj3.eta() - pj2.eta());
      const double alpha = atan(dH/dPhi);
      _histAlpha->fill(alpha*180./PI, weight);
    }


    /// Finalize
    void finalize() {
      /// @todo Apply correction
      // double a, b, c, erra, errb, errc;
      // for (int ibin = 0;  ibin < _histAlpha->getNbins(); ++ibin) {
      // a = _histAlpha->GetBinContent(ibin);
      // erra = _histAlpha->GetBinError(ibin);
      // b = _histAlphaIdeal->GetBinContent(ibin);
      // errb = _histAlphaIdeal->GetBinError(ibin);
      // c = _histAlphaCDF->GetBinContent(ibin);
      // errc = _histAlphaCDF->GetBinError(ibin);
      // _histAlpha->SetBinContent(ibin, b/c);
      // _histAlpha->SetBinError(ibin, sqrt(sqr(b)/sqr(c)*sqr(erra) + sqr(a)/sqr(c)*sqr(errb) +
      // sqr(a*b/(sqr(c)))*sqr(errc) ) );
      // }
      /// @todo Same correction to be applied for _hisR23 and _histJet3eta histograms

      //getLog() << Log::INFO << "Cross-section = " << crossSection()/picobarn << " pb" << endl;
      normalize(_histJet1Et);
      normalize(_histJet2Et);
      normalize(_histR23);
      normalize(_histJet3eta);
      normalize(_histAlpha);
    }

    //@}


  private:

    /// @name Histogram collections
    //@{
    // AIDA::IHistogram2D* _histHvsDphi;
    // AIDA::IHistogram2D* _histRvsAlpha;
    AIDA::IHistogram1D* _histJet1Et;
    AIDA::IHistogram1D* _histJet2Et;
    AIDA::IHistogram1D* _histR23;
    AIDA::IHistogram1D* _histJet3eta;
    AIDA::IHistogram1D* _histAlpha;
    // AIDA::IHistogram1D* _histAlphaMCvsDat;
    AIDA::IHistogram1D* _histAlphaIdeal;
    AIDA::IHistogram1D* _histAlphaCDF;
    AIDA::IHistogram1D* _histR23Ideal;
    AIDA::IHistogram1D* _histR23CDF;
    AIDA::IHistogram1D* _histJet3etaIdeal;
    AIDA::IHistogram1D* _histJet3etaCDF;
    //@}

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CDF_1994_S2952106> plugin_CDF_1994_S2952106;

}
