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
      _histAlphaCDF = bookHistogram1D(6,1,1);

      /// @todo Need better title
      //_histAlphaIdeal = bookHistogram1D("alphaIdeal", 42, -PI/2., PI/2.);
      _histAlphaIdeal = bookHistogram1D(6,1,2);


      /// @todo Need better title
      //_histR23CDF = bookHistogram1D("R23CDF", 50, 0., 5.);
      _histR23CDF = bookHistogram1D(7,1,1);

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

      const double alpha_CDF_sim[] = {0.0517, 0.0461, 0.049, 0.0452, 0.0451, 
				      0.0435, 0.0317, 0.0287, 0.0294, 0.0261, 
				      0.0231, 0.022, 0.0233, 0.0192, 0.0213, 
				      0.0166, 0.0176, 0.0146, 0.0136, 0.0156, 
				      0.0142, 0.0152, 0.0151, 0.0147, 0.0164, 
				      0.0186, 0.018, 0.021, 0.0198, 0.0189, 
				      0.0197, 0.0211, 0.027, 0.0236, 0.0243, 
				      0.0269, 0.0257, 0.0276, 0.0246, 0.0286};
      
      const double alpha_CDF_sim_err[] = {0.0024, 0.0025, 0.0024, 0.0024, 0.0024, 
					 0.0022, 0.0019, 0.0018, 0.0019, 0.0016, 
					 0.0017, 0.0017, 0.0019, 0.0013, 0.0017, 
					 0.0014, 0.0016, 0.0013, 0.0012, 0.0009, 
					 0.0014, 0.0014, 0.0014, 0.0014, 0.0014, 
					 0.0015, 0.0014, 0.0016, 0.0016, 0.0015, 
					 0.0016, 0.0016, 0.0019, 0.0017, 0.0019, 
					 0.0018, 0.0018, 0.0018, 0.0018, 0.0019};

      const double alpha_Ideal_sim[] = {0.0552, 0.0558, 0.0583, 0.055, 0.0495, 
				       0.0433, 0.0393, 0.0346, 0.0331, 0.0296, 
				       0.0258, 0.0196, 0.0171, 0.0179, 0.0174, 
				       0.0141, 0.0114, 0.0096, 0.0076, 0.0087, 
				       0.0099, 0.0079, 0.0102, 0.0114, 0.0124, 
				       0.013, 0.0165, 0.016, 0.0177, 0.019, 
				       0.0232, 0.0243, 0.0238, 0.0248, 0.0235, 
				       0.0298, 0.0292, 0.0291, 0.0268, 0.0316};

      for (int ibin = 0;  ibin < 40; ++ibin) {
	_histAlphaCDF->fill(alpha_CDF_sim[ibin], 
			    alpha_CDF_sim_err[ibin]*alpha_CDF_sim_err[ibin]/alpha_CDF_sim[ibin]);

      }

      /*
      for (int ibin = 0;  ibin < 40; ++ibin) {
	double raw = _histAlpha->GetBinContent(ibin);
	double corr = alpha_CDF_sim[ibin]/alpha_Ideal_sim[ibin];
	_histAlpha->SetBinContent(raw*corr);
	double rawerr = _histAlpha->GetBinError(ibin);
	double CDF_sim_err = alpha_CDF_sim_err[ibin];
	_histAlpha->SetBinError(ibin, sqrt(corr*corr*rawerr*rawerr + 
					   raw*raw/(alpha_Ideal_sim[ibin]*alpha_Ideal_sim[ibin])*
					   alpha_CDF_sim_err[ibin]*alpha_CDF_sim_err[ibin]) );
      }
      */

      const double R23_CDF_sim[] = {0.0005, 0.0161, 0.057, 0.0762, 0.0723, 
				   0.0705, 0.0598, 0.0563, 0.0557, 0.0579, 
				   0.0538, 0.0522, 0.0486, 0.0449, 0.0418, 
				   0.0361, 0.0326, 0.0304, 0.0252, 0.0212, 
				   0.0173, 0.0176, 0.0145, 0.0127, 0.0103, 
				   0.0065, 0.0049, 0.0045, 0.0035, 0.0029, 
				   0.0024, 0.0014, 0.0011, 0.001, 0.0009};

      const double R23_CDF_sim_err[] = {0.0013, 0.0009, 0.0022, 0.0029, 0.0026, 
				       0.0024, 0.0022, 0.0025, 0.0023, 0.0024, 
				       0.0021, 0.0021, 0.0021, 0.0021, 0.0021, 
				       0.0019, 0.0019, 0.0016, 0.0017, 0.0014, 
				       0.001, 0.0014, 0.0012, 0.0013, 0.001, 
				       0.0011, 0.001, 0.001, 0.001, 0.0011, 
				       0.0011, 0.0009, 0.0008, 0.0008, 0.0009};

      const double R23_Ideal_sim[] = {0.0005, 0.0176, 0.0585, 0.0862, 0.0843, 
				     0.0756, 0.0673, 0.0635, 0.0586, 0.0619, 
				     0.0565, 0.0515, 0.0466, 0.0472, 0.0349, 
				     0.0349, 0.0266, 0.0254, 0.0204, 0.0179, 
				     0.0142, 0.0134, 0.0101, 0.009, 0.008, 
				     0.0034, 0.003, 0.0033, 0.0027, 0.0021, 
				     0.0012, 0.0006, 0.0004, 0.0005, 0.0003};
      
      const double eta3_CDF_sim[] = {0.0013, 0.0037, 0.0047, 0.0071, 0.0093, 
				    0.0117, 0.0151, 0.0149, 0.0197, 0.0257, 
				    0.0344, 0.0409, 0.0481, 0.0454, 0.0394, 
				    0.0409, 0.0387, 0.0387, 0.0322, 0.0313, 
				    0.029, 0.0309, 0.0412, 0.0417, 0.0412, 
				    0.0397, 0.0417, 0.0414, 0.0376, 0.0316, 
				    0.027, 0.0186, 0.0186, 0.0132, 0.0127, 
				    0.0106, 0.0071, 0.004, 0.002, 0.0013};

      const double eta3_CDF_sim_err[] = {0.0009, 0.0009, 0.0007, 0.0007, 0.0007, 
					0.001, 0.0012, 0.0012, 0.0013, 0.0016, 
					0.0017, 0.002, 0.002, 0.0022, 0.002, 
					0.002, 0.0018, 0.0018, 0.0016, 0.0017, 
					0.0017, 0.0019, 0.002, 0.0021, 0.002, 
					0.002, 0.0019, 0.002, 0.0018, 0.0017, 
					0.0017, 0.0014, 0.0014, 0.0009, 0.001, 
					0.0009, 0.0009, 0.0008, 0.0008, 0.0009};

      const double eta3_Ideal_sim[] = {0.0017, 0.003, 0.0033, 0.0062, 0.0062, 
				      0.0112, 0.0177, 0.0164, 0.0196, 0.0274, 
				      0.0351, 0.0413, 0.052, 0.0497, 0.0448, 
				      0.0446, 0.0375, 0.0329, 0.0291, 0.0272, 
				      0.0233, 0.0288, 0.0384, 0.0396, 0.0468, 
				      0.0419, 0.0459, 0.0399, 0.0355, 0.0329, 
				      0.0274, 0.023, 0.0201, 0.012, 0.01, 
				      0.008, 0.0051, 0.0051, 0.001, 0.001};




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
      normalize(_histJet1Et, 12.3025); //norm to integral of Ref data
      normalize(_histJet2Et, 12.4565); //norm to integral of Ref data
      normalize(_histJet3eta, 0.19864); //norm to integral of Ref data
      normalize(_histR23, 0.125675); //norm to integral of Ref data
      normalize(_histAlpha, 4.5099); //norm to integral of Ref data
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
