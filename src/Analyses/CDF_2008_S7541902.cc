// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2008_S7541902.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/SVertex.hh"
#include <algorithm>

namespace Rivet {


  // Constructor
  CDF_2008_S7541902::CDF_2008_S7541902()
    : _electronETCut(20.0*GeV), _electronETACut(1.1),
      _eTmissCut(30.0*GeV), _mTCut(20.0*GeV),
      _jetEtCutA(20.0*GeV),  _jetEtCutB(25.0*GeV), _jetETA(2.0),
      _xpoint(1960.)
  {
    setBeams(PROTON, ANTIPROTON);
    setNeedsCrossSection(true);
    
    // Basic FS
    FinalState fs(-3.6, 3.6);
    addProjection(fs, "FS");
    
    // Create a final state with any e-nu pair with invariant mass 65 -> 95 GeV and ET > 20 (W decay products)
    std::vector<std::pair<long,long> > vids;
    vids.push_back(make_pair(ELECTRON, NU_EBAR));
    vids.push_back(make_pair(POSITRON, NU_E));
    FinalState fs2(-3.6, 3.6, 20*GeV);
    InvMassFinalState invfs(fs2, vids, 65*GeV, 95*GeV);
    addProjection(invfs, "INVFS");
    
    // Make a final state without the W decay products for jet clustering
    VetoedFinalState vfs(fs);
    vfs.addVetoOnThisFinalState(invfs);
    addProjection(vfs, "VFS");
    addProjection(FastJets(vfs, FastJets::CDFJETCLU, 0.4), "Jets");
  }


  void CDF_2008_S7541902::init() {
     for (int i = 0 ; i < 4 ; ++i) {
      stringstream title;
      title << "$E_\\perp$ of jet #" << i+1;
      _histJetEt[i] = bookHistogram1D(i+1, 1, 1, title.str(),"$E_\\perp$ [GeV]","$d\\sigma/dE_{T} [pb/GeV]$");
      stringstream title2;
      title2 << "$\\sigma(" << i+1 << " \\text{ jets})/\\sigma(" << i <<" \\text{ jets})$" ;
      _histJetMultRatio[i] = bookDataPointSet(5 , 1, i+1, title2.str(),"$\\sqrt(s)$ [GeV]",title2.str());
      stringstream title3;
      title3 << "$\\sigma(" << i+1 << " \\text{ jets})$" ;
      _histJetMult[i]   = bookHistogram1D(i+6, 1, 1, title3.str(),"$\\sqrt(s)$ [GeV]",title3.str());
     } 
     _histJetMultNorm = bookHistogram1D("norm", "$\\sigma(0 \\text{ jets})$", "$\\sqrt(s)$ [GeV]", "No. of events", 1, _xpoint, _xpoint+1.);
  }



  void CDF_2008_S7541902::analyze(const Event& event) {
    // Get the W decay products (electron and neutrino)
    const InvMassFinalState& invMassFinalState = applyProjection<InvMassFinalState>(event, "INVFS");
    const ParticleVector&  wDecayProducts = invMassFinalState.particles();
        
    FourMomentum electronP, neutrinoP;
    bool gotElectron(false), gotNeutrino(false);
    foreach (const Particle& p, wDecayProducts) {
      FourMomentum p4 = p.momentum();
      if (p4.Et() > _electronETCut && fabs(p4.eta()) < _electronETACut && abs(p.pdgId()) == ELECTRON) {
        electronP = p4;
        gotElectron = true;
      }
      else if (p4.Et() > _eTmissCut && abs(p.pdgId()) == NU_E) {
        neutrinoP = p4;
        gotNeutrino = true;
      }
    }
    
    // Veto event if the electron or MET cuts fail
    if (!gotElectron || !gotNeutrino) vetoEvent;

    // Veto event if the MTR cut fails
    double mT2 = 2.0 * ( electronP.pT()*neutrinoP.pT() - electronP.px()*neutrinoP.px() - electronP.py()*neutrinoP.py() );
    if (sqrt(mT2) < _mTCut ) vetoEvent;

    // Get the jets
    const JetAlg& jetProj = applyProjection<FastJets>(event, "Jets");
    Jets theJets = jetProj.jetsByEt(_jetEtCutA);
    size_t njetsA(0), njetsB(0);
    foreach (const Jet& j, theJets) {
      const FourMomentum pj = j.momentum();
      if (fabs(pj.rapidity()) < _jetETA) {
        // Fill differential histograms for top 4 jets with Et > 20
        if (njetsA < 4 && pj.Et() > _jetEtCutA) {
          ++njetsA;
          _histJetEt[njetsA-1]->fill(pj.Et(), event.weight());
        }
        // Count number of jets with Et > 25 (for multiplicity histograms)
        if (pj.Et() > _jetEtCutB) ++njetsB;
      }
    }
 
    // Jet multiplicity
    _histJetMultNorm->fill(_xpoint, event.weight());
    for (size_t i = 1; i <= njetsB ; ++i) {
      _histJetMult[i-1]->fill(_xpoint, event.weight());
    }
  }
  
  
  
  // Finalize
  void CDF_2008_S7541902::finalize() { 
    float xsec = crossSection()/sumOfWeights();
    // Get the x-axis for the ratio plots
    /// @todo Replace with autobooking etc. once YODA in place    
    std::vector<double> xval; xval.push_back(_xpoint);
    std::vector<double> xerr; xerr.push_back(.5);
    // Fill the first ratio histogram using the special normalisation histogram for the total cross section
    double ratio1to0 = 0.;
    if (_histJetMultNorm->binHeight(0) > 0.) ratio1to0 = _histJetMult[0]->binHeight(0)/_histJetMultNorm->binHeight(0);
    // Get the fractional error on the ratio histogram
    double frac_err1to0 = 0.;
    if (_histJetMult[0]->binHeight(0) > 0.)  frac_err1to0 = _histJetMult[0]->binError(0)/_histJetMult[0]->binHeight(0);
    if (_histJetMultNorm->binHeight(0) > 0.) {
      frac_err1to0 *= frac_err1to0;
      frac_err1to0 += pow(_histJetMultNorm->binError(0)/_histJetMultNorm->binHeight(0),2.);
      frac_err1to0 = sqrt(frac_err1to0);
    }

    /// @todo Replace with autobooking etc. once YODA in place    
    std::vector<double> yval[4]; yval[0].push_back(ratio1to0);
    std::vector<double> yerr[4]; yerr[0].push_back(ratio1to0*frac_err1to0);
    _histJetMultRatio[0]->setCoordinate(0,xval,xerr);
    _histJetMultRatio[0]->setCoordinate(1,yval[0],yerr[0]);
    for (int i = 0; i < 4; ++i) {
      if (i < 3) {
        float ratio = 0.0;
        if (_histJetMult[i]->binHeight(0) > 0.0) ratio = _histJetMult[i+1]->binHeight(0)/_histJetMult[i]->binHeight(0);
        float frac_err = 0.0;
        if (_histJetMult[i]->binHeight(0) > 0.0) frac_err = _histJetMult[i]->binError(0)/_histJetMult[i]->binHeight(0);
        if (_histJetMult[i+1]->binHeight(0) > 0.0) {
          frac_err *= frac_err;
          frac_err += pow(_histJetMult[i+1]->binError(0)/_histJetMult[i+1]->binHeight(0),2.);
          frac_err = sqrt(frac_err);
        }
        yval[i+1].push_back(ratio);
        yerr[i+1].push_back(ratio*frac_err);
        _histJetMultRatio[i+1]->setCoordinate(0,xval,xerr);
        _histJetMultRatio[i+1]->setCoordinate(1,yval[i+1],yerr[i+1]);
      }
      _histJetEt[i]->scale(xsec);
      _histJetMult[i]->scale(xsec);
    }
    _histJetMultNorm->scale(xsec);
  }
}
