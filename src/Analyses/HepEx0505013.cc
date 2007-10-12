// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/HepEx0505013.hh"
#include "Rivet/RivetAIDA.hh"


namespace Rivet {
  
  // Book histograms
  void HepEx0505013::init() {

    string hist_title[18][2];
    
    for (int i=0; i<6; ++i) { //18=6x3 pT bins, one histogram each
      for (int j=0; j<3; ++j) {
	int k = i*3+j;
	stringstream lStream; 
	lStream << k+1;
	hist_title[k][0] = "Differential Jet Shape Rho, Pt bin " + lStream.str();
	_histRho_pT[k]  = bookHistogram1D(i+1, 1, j+1, hist_title[k][0]);
	
	hist_title[k][1] = "Integral Jet Shape Psi, Pt bin " + lStream.str();
	_histPsi_pT[k]  = bookHistogram1D(6+i+1, 1, j+1, hist_title[k][1]);
      }
    }
    
    string hist_title_oneminPsi = "One minus Psi(0.3 over R)";
    _histOneMinPsi = bookHistogram1D(13, 1, 1, hist_title_oneminPsi);
    
  }
  
  
  // Do the analysis
  void HepEx0505013::analyze(const Event& event) {
    Log log = getLog();
    log << Log::DEBUG << "Starting analyzing" << endl;
    
    // Analyse and print some info  
    #ifdef HAVE_FASTJET
    const FastJets& jetpro = event.applyProjection(_jetsproj);
    #else
    const D0ILConeJets& jetpro = event.applyProjection(_jetsproj);
    ////const KtJets& jetpro = event.applyProjection(_jetsproj);    
    #endif
    
    log << Log::DEBUG << "Jet multiplicity before any pT cut = " << jetpro.getNJets() << endl;
    
    
    // Find vertex and check  that its z-component is < 60 cm from the nominal IP
    //const PVertex& pv = event.applyProjection(_vertexproj);
    /// @todo SEGV: either the HepMC event record is not filled properly or the F77-Wrapper functions are faulty
    /// @todo z- value assumed to be in mm, PYTHIA convention: dangerous!
    //if (fabs(pv.getPrimaryVertex().position().z()) < 600.0) {


    /// @todo Remove compile-time flags like this when possible.
    #ifdef HAVE_FASTJET
    typedef vector<fastjet::PseudoJet> Jets;
    const Jets& jets = jetpro.getJetsPt();
    #else
    typedef list<LorentzVector> Jets;
    const Jets& jets = jetpro.getLorentzJets();
    ////const vector<KtJet::KtLorentzVector>& jets = jetpro.getJets();
    #endif
    
    log << Log::DEBUG << "jetlist size = " << jets.size() << endl;
    
    int Njet = 0;
    bool jetcutpass = false;
    
    for (Jets::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {
      log << Log::DEBUG << "List item pT = " << jt->perp() << " E=" << jt->e() << " pz=" << jt->pz() << endl;
      if (jt->perp()>37. && fabs(jt->rapidity())>0.1 && fabs(jt->rapidity())<0.7) jetcutpass = true;
      ++Njet;
      log << Log::DEBUG << "Jet pT =" << jt->perp() << " y=" << jt->rapidity() << " phi=" << jt->phi() << endl; 
    }
    
    if (jetcutpass) {
      
      const TotalVisibleMomentum& caloMissEt = event.applyProjection(_calmetproj);
      log << Log::DEBUG << "CaloMissEt.getMomentum().perp() = " << caloMissEt.getMomentum().perp() << endl;

      if (caloMissEt.getMomentum().perp()/sqrt(caloMissEt.getSET()) < 3.5) {
        
        LorentzVector jetaxis;
        _jetaxes.clear();
        for (Jets::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {

          if (fabs(jt->rapidity())<1.1) { //only Central Calorimeter jets
            jetaxis.setPx(jt->px());
            jetaxis.setPy(jt->py());
            jetaxis.setPz(jt->pz());
            jetaxis.setE(jt->e());

            _jetaxes.push_back(jetaxis);
          }
        }
        if (_jetaxes.size()>0) { //determine jet shapes
          
	  
          const JetShape& jetShape = event.applyProjection(_jetshapeproj);
	  
          //fill histograms
          for (unsigned int jind=0; jind<_jetaxes.size(); ++jind) {
            for (int ipT=0; ipT<18; ++ipT) {
              if (_jetaxes[jind].perp() > _pTbins[ipT] && _jetaxes[jind].perp() <= _pTbins[ipT+1]) {
                _ShapeWeights[ipT] += event.weight(); 
                for (int rbin=0; rbin<_jetshapeproj.getNbins(); ++rbin) {
                  double rad = _jetshapeproj.getRmin() +(rbin+0.5)*_jetshapeproj.getInterval();
                  _histRho_pT[ipT]->fill(rad/_Rjet, double(_diffjetshapes[jind][rbin]) * event.weight());
                  _histPsi_pT[ipT]->fill(rad/_Rjet, double(_intjetshapes[jind][rbin]) * event.weight());
                }
                _histOneMinPsi->fill((_pTbins[ipT]+_pTbins[ipT+1])/2., 
                                     _oneminPsiShape[jind] * event.weight());
              }
            }
          }	  
        }
      }
    }
    
    
    
    // Finished
    log << Log::DEBUG << "Finished analyzing" << endl;
  }
  
  
  // Finalize
  void HepEx0505013::finalize() { 
    
    vector<double> OneMinPsiBins(18, 0.);
    
    for (int i=0; i<18; ++i) { //18 pT bins, one histogram each
      
      if (_ShapeWeights[i] > 0.) {
        //divide each histogram, each bin, by sum of event weights _EventWeights
        _histRho_pT[i]->scale(1./_ShapeWeights[i]);
        _histPsi_pT[i]->scale(1./_ShapeWeights[i]);	
        OneMinPsiBins[i] = _histOneMinPsi->binHeight(i)/_ShapeWeights[i];
      }
      else OneMinPsiBins[i] = 0.; 
    }
    
    _histOneMinPsi->reset();
    for (int i=0; i<18; ++i) { //refill the reset histogram with normalized bins
      _histOneMinPsi->fill((_pTbins[i]+_pTbins[i+1])/2., OneMinPsiBins[i]);
    }
    
  }
  
  
}
