// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/CDF_1994_S2952106.hh"

namespace Rivet {


  void CDF_1994_S2952106::init() {
    /// @todo Use histogram auto-booking

    //const string hname = "HvsDphi";
    //const string htitle = "H vs Delta phi";
    //_histHvsDphi = bookHistogram2D(hname, htitle, 40, -4., 4., 32, 0., 3.2);

    //const string hname2 = "RvsAlpha";
    //const string htitle2 = "R vs alpha";
    //_histRvsAlpha = bookHistogram2D(hname2, htitle2, 50, 0., 5., 32, -1.6, 1.6);

    const string hname3 = "Jet1Et";
    const string htitle3 = "$E_\\perp$ of leading jet";
    _histJet1Et = bookHistogram1D(hname3, htitle3, 40, 0., 500.);

    const string hname4 = "Jet2Et";
    const string htitle4 = "$E_\\perp$ of 2nd leading jet";
    _histJet2Et = bookHistogram1D(hname4, htitle4, 40, 0., 500.);

    const string hname5 = "R23";
    const string htitle5 = "$R$ distance between 2nd and 3rd jet";
    _histR23 = bookHistogram1D(hname5, htitle5, 50, 0., 5.);

    const string hname6 = "Jet3eta";
    const string htitle6 = "Pseudorapidity, $\\eta$, of 3rd jet";
    _histJet3eta = bookHistogram1D(hname6, htitle6, 42, -4., 4.);

    /// @todo Need better title
    const string hname7 = "alpha";
    const string htitle7 = "$\\alpha$";
    _histAlpha = bookHistogram1D(hname7, htitle7, 42, -PI/2., PI/2.);

    //const string hname8 = "alphaMCvsDat";
    //const string htitle8 = "alpha MC vs. Data ";
    //_histAlphaMCvsDat = bookHistogram1D(hname8, htitle8, 42, -PI/2., PI/2.);

    /// @todo Need better title
    const string hname9 = "alphaIdeal";
    const string htitle9 = "$\\alpha_\\text{ideal}$";
    _histAlpaIdeal = bookHistogram1D(hname9, htitle9, 42, -PI/2., PI/2.);

    /// @todo Need better title
    const string hname10 = "alphaCDF";
    const string htitle10 = "$\\alpha_\\text{CDF}$";
    _histAlpaCDF = bookHistogram1D(hname10, htitle10, 42, -PI/2., PI/2.);

    /// @todo Need better title
    const string hname11 = "R23Ideal";
    const string htitle11 = "$R_{23}^\\text{ideal}$";  
    _histR23Ideal = bookHistogram1D(hname11, htitle11, 50, 0., 5.);

    /// @todo Need better title
    const string hname12 = "R23CDF";
    const string htitle12 = "$R_{23}^\\text{CDF}$";
    _histR23CDF = bookHistogram1D(hname12, htitle12, 50, 0., 5.);

    /// @todo Need better title
    const string hname13 = "Jet3etaIdeal";
    const string htitle13 = "Jet #3 $\\eta_\\text{ideal}$";  
    _histJet3etaIdeal = bookHistogram1D(hname13, htitle13, 42, -4., 4.);

    /// @todo Need better title
    const string hname14 = "Jet3etaCDF";
    const string htitle14 = "Jet #3 $\\eta_\\text{CDF}$";  
    _histJet3etaCDF = bookHistogram1D(hname14, htitle14, 42, -4., 4.);

    /// @todo Move to constructor
    _events3jPassed = 0.0;
  }



  // Do the analysis
  void CDF_1994_S2952106::analyze(const Event & event) {
    const FastJets& jetpro = applyProjection<FastJets>(event, "ConeJets");
    getLog() << Log::DEBUG << "Jet multiplicity before any pT cut = " << jetpro.size() << endl;

    _fail = true;

    // Find vertex and check  that its z-component is < 60 cm from the nominal IP
    const PVertex& pv = applyProjection<PVertex>(event, "PV");
    if (fabs(pv.getPVPosition().z())/mm < _pvzmax) {
      const TotalVisibleMomentum& caloMissEt = applyProjection<TotalVisibleMomentum>(event, "CalMET");
      getLog() << Log::DEBUG << "CaloMissEt.getMomentum().pT() = " << caloMissEt.getMomentum().pT() << endl;
      if (caloMissEt.getMomentum().pT()/sqrt(caloMissEt.getSET()) < _metsetmax) {
        PseudoJets jets = jetpro.pseudoJets();
        PseudoJets::const_iterator jet1stPt = jets.end();
        PseudoJets::const_iterator jet2ndPt = jets.end();
        PseudoJets::const_iterator jet3rdPt = jets.end();
        getLog() << Log::DEBUG << "jetlist size = " << jets.size() << endl;

        int Njet = 0;
        int NjetPtCut = 0;
        for (PseudoJets::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {
          getLog() << Log::DEBUG << "List item pT = " << jt->perp() << " E=" << jt->E() 
              << " pz=" << jt->pz() << endl;
          if (jt->perp() > _leadJetPt) ++NjetPtCut;
          getLog() << Log::DEBUG << "Jet pT =" << jt->perp() << " y=" << jt->pseudorapidity() 
              << " phi=" << jt->phi() << endl; 
          if (jet1stPt == jets.end() || jt->perp() > jet1stPt->perp()) {
            jet3rdPt = jet2ndPt;
            jet2ndPt = jet1stPt;
            jet1stPt = jt;
            Njet++;
          } else if (jet2ndPt == jets.end() || jt->perp() > jet2ndPt->perp()) {
            jet3rdPt = jet2ndPt;
            jet2ndPt = jt;
            Njet++;
          } else if (jet3rdPt == jets.end() || jt->perp() > jet3rdPt->perp()) {
            jet3rdPt = jt;
            Njet++;
          }
        }
        
        if (NjetPtCut >= 1) {
          getLog() << Log::DEBUG << "Jet multiplicity after pT > 100GeV cut = " << NjetPtCut << endl; 
        }
        
        if (Njet>=3 && fabs(jet1stPt->pseudorapidity())<_etamax && fabs(jet2ndPt->pseudorapidity())<_etamax) {
          getLog() << Log::DEBUG << "Jet eta and pT requirements fulfilled" << endl;
          
          if (fabs(fabs(jet1stPt->phi()-jet2ndPt->phi())-PI) < _phimin) {
            getLog() << Log::DEBUG << "1st & 2nd Jet phi requirement fulfilled" << endl;
            
            _fail = false;

            _histJet1Et->fill(jet1stPt->perp(), event.weight());
            _histJet2Et->fill(jet2ndPt->perp(), event.weight());
            //_histR23->fill(jet2ndPt->deltaR(*jet3rdPt), event.weight());
            _histR23->fill(delta_rad(jet2ndPt->pseudorapidity(),jet2ndPt->phi(), 
                                     jet3rdPt->pseudorapidity(),jet3rdPt->phi() ), 
                           event.weight());
            _histJet3eta->fill(jet3rdPt->pseudorapidity(), event.weight());
            
            
            // Next cut only required for alpha studies
            if (jet3rdPt->perp() > _3rdJetPt) {
              getLog() << Log::DEBUG << "jet3rdPt->pT()=" << jet3rdPt->perp() << " (>10.)" << endl;
              
	      _events3jPassed += event.weight();

              double dPhi = fabs(jet3rdPt->phi() - jet2ndPt->phi());
              dPhi -= int(dPhi/PI); //dPhi % PI (modulo)
              
              double dH = jet3rdPt->pseudorapidity() - jet2ndPt->pseudorapidity();
              if (jet2ndPt->pseudorapidity() < 0.) dH *= -1;
              
              double alpha = atan(dH/dPhi);
              //cout << "alpha=" << alpha << "  dH=" << dH 
              // << "  dPhi=" << dPhi << endl;
              
              _histAlpha->fill(alpha, event.weight());
              //_histAlphaMCvsDat->fill(alpha, event.weight());
              
            } //3rd jet Pt cut
          } //delta phi 1st, 2nd jet
        } //1st + 2nd jet pseudoRapidity & Njet>=3
      } //MET/sqrt(SET) cut
    } //z-vertex

    if (_fail) vetoEvent(event); 

  }
  
  
  // Finalize
  void CDF_1994_S2952106::finalize() { 
    /*
   /// @todo Apply correction
   double a, b, c, erra, errb, errc;
   for (int ibin = 0;  ibin < _histAlpha->getNbins(); ++ibin) {
   a = _histAlpha->GetBinContent(ibin);
   erra = _histAlpha->GetBinError(ibin);
   b = _histAlpaIdeal->GetBinContent(ibin);
   errb = _histAlpaIdeal->GetBinError(ibin);
   c = _histAlpaCDF->GetBinContent(ibin);
   errc = _histAlpaCDF->GetBinError(ibin);
   _histAlpha->SetBinContent(ibin, b/c);
   _histAlpha->SetBinError(ibin, sqrt(sqr(b)/sqr(c)*sqr(erra) + sqr(a)/sqr(c)*sqr(errb) + 
   sqr(a*b/(sqr(c)))*sqr(errc) ) );
   }
   /// @todo Same correction to be applied for _hisR23 and _histJet3eta histograms
   */
        
    // Normalise histograms to integrated publication luminosity of 4.2 pb^-1 
    //const double fac = 4.2 / picobarn * crossSection() / sumOfWeights();
    const double fac = 1. / sumOfWeights();
    getLog() << Log::INFO 
             << "Cross-section = " << crossSection()/picobarn << " pb"
             << " -> scale factor = " << fac << endl;
    _histJet1Et->scale(fac);
    _histJet2Et->scale(fac);
    _histR23->scale(fac);
    _histJet3eta->scale(fac);

    const double fac3j = 1. / _events3jPassed;

    _histAlpha->scale(fac3j);
  }
  
  
}
