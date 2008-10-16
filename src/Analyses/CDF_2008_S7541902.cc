// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2008_S7541902.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/AnalysisHandler.hh"
#include <algorithm>

namespace Rivet {


  void CDF_2008_S7541902::init() {

     for (int i = 0 ; i < 4 ; i++) {
      stringstream title;
      title << "Et of jet " << i + 1;
      _histJetEt[i]    = bookHistogram1D(i+1,1,1,title.str().c_str());
      stringstream title2;
      title2 << "sig(" << i+1 << " jets)/sig(" << i <<" jets)" ;
      _histJetMultRatio[i] = bookDataPointSet(5,1,i+1,title2.str().c_str());
      stringstream title3;
      title3 << "sig(" << i+1 << "jets)" ;
      _histJetMult[i]   = bookHistogram1D(i+6,1,1,title3.str().c_str());
     }
     
     _histJetMultNorm = bookHistogram1D("norm","sig(0 jets)",1,_xpoint,_xpoint+1.);
   
  }



  void CDF_2008_S7541902::analyze(const Event& event) {
    Log log = getLog();
    log << Log::DEBUG << "Starting analyzing" << endl;
    // Get the W deday products (electron and neutrino)
    const InvMassFinalState& invMassFinalState = applyProjection<InvMassFinalState>(event, "INVFS");
    const ParticleVector&  WDecayProducts =  invMassFinalState.particles();
        
    FourMomentum electronP;
    bool gotElectron = false;
    FourMomentum neutrinoP;
    bool gotNeutrino = false;
    for(ParticleVector::const_iterator p =  WDecayProducts.begin();
	p !=  WDecayProducts.end(); ++p){
      if(p->momentum().Et() > _electronETCut  && 
	 fabs(p->momentum().pseudorapidity()) < _electronETACut  && 
	 abs(p->getPdgId())==ELECTRON) {
	electronP = p->momentum();
	gotElectron = true;
      }
      if(p->momentum().Et() > _eTmissCut && 
	 abs(p->getPdgId())==NU_E) {
	neutrinoP = p->momentum();
	gotNeutrino = true;
      }
    }

    // Veto event if the electron or MET cuts fail
    if(!gotElectron || !gotNeutrino){
      vetoEvent(event);
      return;
    }
    // Veto event if the MTR cut fails
    double mT2 = electronP.pT() * neutrinoP.pT() - 
    electronP.px() * neutrinoP.px() -
    electronP.py() * neutrinoP.py();
    if(mT2 < _mT2Cut ){
      vetoEvent(event);
      return;
    }

    // get the jets
    const FastJets &jetProj = applyProjection<FastJets>(event, "Jets");
    
    Jets theJets = jetProj.getJetsByPt(_jetEtCutA);
    sort(theJets.begin(), theJets.end(), Particle::byETDescending());
    
    Jets foundJets;
    int njets = 0;    
    for(Jets::const_iterator jIt =theJets.begin(); 
	foundJets.size()!=4 && jIt != theJets.end(); ++jIt){
      if(fabs(jIt->momentum().rapidity()) < _jetETA) {
        // store all jets with ET > 20. for differential histograms 
	if(jIt->momentum().Et() > _jetEtCutA) foundJets.push_back(*jIt);
        // count number of jets with ET > 25. for multiplicity histograms
        if(jIt->momentum().Et() > _jetEtCutB) njets++; 
      }
     }
    // fill jet ET distributions for first four jets
    for(size_t i = 0; i != foundJets.size(); ++i){
      _histJetEt[i]->fill(foundJets[i].momentum().Et(),event.weight());
    }
 
    // jet multiplicity :
    _histJetMultNorm->fill(_xpoint,event.weight());
    for (size_t i = 1; i <= njets ; ++i) { 
      _histJetMult[i-1]->fill(_xpoint,event.weight()); 
    } 
 
   return;
  }
  
  
  
  // Finalize
  void CDF_2008_S7541902::finalize() { 
    float xsec = crossSection()/getHandler().sumOfWeights();
    // get the x-axis for the ratio plots
    std::vector<double> xval; xval.push_back(_xpoint);
    std::vector<double> xerr; xerr.push_back(0.);
    // fill the first ratio histogram using the special normalisation histogram for the total cross section
    double ratio1to0 = 0.;
    if (_histJetMultNorm->binHeight(0) > 0.) ratio1to0 = _histJetMult[0]->binHeight(0)/_histJetMultNorm->binHeight(0);
    // get the fractional error on the ratio histogram
    double frac_err1to0 = 0.;
    if (_histJetMult[0]->binHeight(0) > 0.)  frac_err1to0 = _histJetMult[0]->binError(0)/_histJetMult[0]->binHeight(0);
    if (_histJetMultNorm->binHeight(0) > 0.) {
      frac_err1to0 *= frac_err1to0;
      frac_err1to0 += pow(_histJetMultNorm->binError(0)/_histJetMultNorm->binHeight(0),2.);
      frac_err1to0 = sqrt(frac_err1to0);
    }

    std::vector<double> yval[4]; yval[0].push_back(ratio1to0);
    std::vector<double> yerr[4]; yerr[0].push_back(ratio1to0*frac_err1to0);
    _histJetMultRatio[0]->setCoordinate(0,xval,xerr);
    _histJetMultRatio[0]->setCoordinate(1,yval[0],yerr[0]);
    for (int i=0; i<4; i++) {
      if(i<3) {
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
