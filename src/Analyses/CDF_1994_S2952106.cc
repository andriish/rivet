// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/CDF_1994_S2952106.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PVertex.hh"

namespace Rivet {


  CDF_1994_S2952106::CDF_1994_S2952106()
    : Analysis("CDF_1994_S2952106"), 
      _pvzmax(600*mm), _leadJetPt(100*GeV), _3rdJetPt(10*GeV),
      _etamax(0.7), _phimin(PI/18.0), _metsetmax(6.0*GeV)
  {
    setBeams(PROTON, ANTIPROTON);
    setNeedsCrossSection(true);

    const FinalState fs(-4.2, 4.2);
    addProjection(fs, "FS");
    addProjection(FastJets(fs, FastJets::CDFJETCLU, 0.7), "ConeJets");
    addProjection(TotalVisibleMomentum(fs), "CalMET");
    addProjection(PVertex(), "PV");

    // Veto (anti)neutrinos, and muons with pT above 1.0 GeV
    VetoedFinalState vfs(fs);
    vfs
      .addVetoPairId(NU_E)
      .addVetoPairId(NU_MU)
      .addVetoPairId(NU_TAU)
      .addVetoDetail(MUON, 1.0*GeV, MAXDOUBLE);
    addProjection(vfs, "VFS");

    _events3jPassed = 0.0;
  }


  void CDF_1994_S2952106::init() {
    /// @todo Use histogram auto-booking

    //const string hname = "HvsDphi";
    //const string htitle = "H vs Delta phi";
    //_histHvsDphi = bookHistogram2D(hname, htitle, 40, -4., 4., 32, 0., 3.2);

    //const string hname2 = "RvsAlpha";
    //const string htitle2 = "R vs alpha";
    //_histRvsAlpha = bookHistogram2D(hname2, htitle2, 50, 0., 5., 32, -1.6, 1.6);

    const string ylabel = "Fraction of events";

    const string hname3 = "Jet1Et";
    const string htitle3 = "$E_\\perp$ of leading jet";
    _histJet1Et = 
      bookHistogram1D(hname3, htitle3, "$E_\\perp^1$ [GeV]", ylabel, 40, 0., 500.);

    const string hname4 = "Jet2Et";
    const string htitle4 = "$E_\\perp$ of 2nd leading jet";
    _histJet2Et = 
      bookHistogram1D(hname4, htitle4, "$E_\\perp^2$ [GeV]", ylabel, 40, 0., 500.);

    const string hname5 = "R23";
    const string htitle5 = "$R$ distance between 2nd and 3rd jet";
    _histR23 = 
      bookHistogram1D(hname5, htitle5, "$R_{23}$", ylabel, 50, 0., 5.);

    const string hname6 = "Jet3eta";
    const string htitle6 = "Pseudorapidity, $\\eta$, of 3rd jet";
    _histJet3eta =
      bookHistogram1D(hname6, htitle6, "$\\eta_3$", ylabel, 42, -4., 4.);

    /// @todo Need better title
    const string hname7 = "alpha";
    const string htitle7 = "$\\alpha$";
    _histAlpha = 
      bookHistogram1D(hname7, htitle7, "$\\alpha$", ylabel, 42, -PI/2., PI/2.);

    //const string hname8 = "alphaMCvsDat";
    //const string htitle8 = "alpha MC vs. Data ";
    //_histAlphaMCvsDat = bookHistogram1D(hname8, htitle8, 42, -PI/2., PI/2.);

    /// @todo Need better title
    const string hname9 = "alphaIdeal";
    const string htitle9 = "$\\alpha_\\text{ideal}$";
    _histAlpaIdeal = 
      bookHistogram1D(hname9, htitle9, "$\\alpha$", ylabel, 42, -PI/2., PI/2.);

    /// @todo Need better title
    const string hname10 = "alphaCDF";
    const string htitle10 = "$\\alpha_\\text{CDF}$";
    _histAlpaCDF = 
      bookHistogram1D(hname10, htitle10, "$\\alpha$", ylabel, 42, -PI/2., PI/2.);

    /// @todo Need better title
    const string hname11 = "R23Ideal";
    const string htitle11 = "$R_{23}^\\text{ideal}$";  
    _histR23Ideal = 
      bookHistogram1D(hname11, htitle11, "$R_{23}$", ylabel, 50, 0., 5.);

    /// @todo Need better title
    const string hname12 = "R23CDF";
    const string htitle12 = "$R_{23}^\\text{CDF}$";
    _histR23CDF = 
      bookHistogram1D(hname12, htitle12, "$R_{23}$", ylabel, 50, 0., 5.);

    /// @todo Need better title
    const string hname13 = "Jet3etaIdeal";
    const string htitle13 = "Jet #3 $\\eta_\\text{ideal}$";  
    _histJet3etaIdeal = 
      bookHistogram1D(hname13, htitle13, "$\\eta_3$", ylabel, 42, -4., 4.);

    /// @todo Need better title
    const string hname14 = "Jet3etaCDF";
    const string htitle14 = "Jet #3 $\\eta_\\text{CDF}$";  
    _histJet3etaCDF = 
      bookHistogram1D(hname14, htitle14, "$\\eta_3$", ylabel, 42, -4., 4.);
  }



  // Do the analysis
  void CDF_1994_S2952106::analyze(const Event & event) {
    const Jets jets = applyProjection<FastJets>(event, "ConeJets").jetsByPt();
    getLog() << Log::DEBUG << "Jet multiplicity before any cuts = " << jets.size() << endl;

    // Find vertex and check  that its z-component is < 60 cm from the nominal IP
    const PVertex& pv = applyProjection<PVertex>(event, "PV");
    if (fabs(pv.position().z())/mm > _pvzmax) {
      vetoEvent;
    }

    // Check there isn't too much missing Et
    const TotalVisibleMomentum& caloMissEt = applyProjection<TotalVisibleMomentum>(event, "CalMET");
    getLog() << Log::DEBUG << "Missing pT = " << caloMissEt.momentum().pT()/GeV << " GeV" << endl;      
    if ((caloMissEt.momentum().pT()/GeV) / sqrt(caloMissEt.scalarET()/GeV) > _metsetmax) {
      vetoEvent;
    }

    // Check jet requirements
    if (jets.size() < 3) vetoEvent;
    if (jets[0].momentum().pT()/GeV < 100) vetoEvent;

    // More jet 1,2,3 checks
    FourMomentum pj1(jets[0].momentum()), pj2(jets[1].momentum()), pj3(jets[2].momentum());
    if (fabs(pj1.eta()) > _etamax || fabs(pj2.eta()) > _etamax) vetoEvent;
    getLog() << Log::DEBUG << "Jet 1 & 2 eta, pT requirements fulfilled" << endl;          

    if (deltaPhi(pj1.phi(), pj2.phi()) > _phimin) vetoEvent;
    getLog() << Log::DEBUG << "Jet 1 & 2 phi requirement fulfilled" << endl;

    const double weight = event.weight();
    _histJet1Et->fill(pj1.pT(), weight);
    _histJet2Et->fill(pj2.pT(), weight);
    _histR23->fill(deltaR(pj2, pj3), weight);
    _histJet3eta->fill(pj3.eta(), weight);
                        
    // Next cut only required for alpha studies
    if (pj3.pT() < _3rdJetPt) vetoEvent;
    getLog() << Log::DEBUG << "3rd jet passes alpha histo pT cut" << endl;      
    _events3jPassed += weight;

    // Calc and plot alpha
    const double dPhi = deltaPhi(pj3.phi(), pj2.phi());    
    const double dH = sign(pj2.eta()) * (pj3.eta() - pj2.eta());
    const double alpha = atan(dH/dPhi);
    _histAlpha->fill(alpha, weight);
  }
  
  
  // Finalize
  void CDF_1994_S2952106::finalize() { 
   /// @todo Apply correction
   // double a, b, c, erra, errb, errc;
   // for (int ibin = 0;  ibin < _histAlpha->getNbins(); ++ibin) {
   // a = _histAlpha->GetBinContent(ibin);
   // erra = _histAlpha->GetBinError(ibin);
   // b = _histAlpaIdeal->GetBinContent(ibin);
   // errb = _histAlpaIdeal->GetBinError(ibin);
   // c = _histAlpaCDF->GetBinContent(ibin);
   // errc = _histAlpaCDF->GetBinError(ibin);
   // _histAlpha->SetBinContent(ibin, b/c);
   // _histAlpha->SetBinError(ibin, sqrt(sqr(b)/sqr(c)*sqr(erra) + sqr(a)/sqr(c)*sqr(errb) + 
   // sqr(a*b/(sqr(c)))*sqr(errc) ) );
   // }
   /// @todo Same correction to be applied for _hisR23 and _histJet3eta histograms
        
    getLog() << Log::INFO << "Cross-section = " << crossSection()/picobarn << " pb" << endl;
    normalize(_histJet1Et);
    normalize(_histJet2Et);
    normalize(_histR23);
    normalize(_histJet3eta);
    normalize(_histAlpha);
  }
  
  
}
