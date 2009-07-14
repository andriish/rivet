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

    _histJet1Et  = bookHistogram1D("Jet1Et", 40, 0., 500.);
    _histJet2Et  = bookHistogram1D("Jet2Et", 40, 0., 500.);
    _histR23     = bookHistogram1D("R23", 50, 0., 5.);
    _histJet3eta = bookHistogram1D("Jet3eta", 42, -4., 4.);

    /// @todo Need better title
    _histAlpha = bookHistogram1D("alpha", 42, -PI/2., PI/2.);

    //const string hname8 = "alphaMCvsDat";
    //const string htitle8 = "alpha MC vs. Data ";
    //_histAlphaMCvsDat = bookHistogram1D(hname8, htitle8, 42, -PI/2., PI/2.);

    /// @todo Need better title
    _histAlpaIdeal = bookHistogram1D("alphaIdeal", 42, -PI/2., PI/2.);

    /// @todo Need better title
    _histAlpaCDF = bookHistogram1D("alphaCDF", 42, -PI/2., PI/2.);

    /// @todo Need better title
    _histR23Ideal = bookHistogram1D("R23Ideal", 50, 0., 5.);

    /// @todo Need better title
    _histR23CDF = bookHistogram1D("R23CDF", 50, 0., 5.);

    /// @todo Need better title
    _histJet3etaIdeal = bookHistogram1D("Jet3etaIdeal", 42, -4., 4.);

    /// @todo Need better title
    _histJet3etaCDF = bookHistogram1D("Jet3etaCDF", 42, -4., 4.);
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
