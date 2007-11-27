// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2005_S6217184.hh"
#include "Rivet/RivetAIDA.hh"


namespace Rivet {
  
  // Book histograms
  void CDF_2005_S6217184::init() {

    string hist_title[18][2];
    
    for (int i=0; i<6; ++i) { //18=6x3 pT bins, one histogram each
      for (int j=0; j<3; ++j) {
        int k = i*3+j;
        stringstream lStream; 
        lStream << k+1;
        hist_title[k][0] = "Differential Jet Shape Rho, Pt bin " + lStream.str();
        _profhistRho_pT[k]  = bookProfile1D(i+1, 1, j+1, hist_title[k][0]);
        
        hist_title[k][1] = "Integral Jet Shape Psi, Pt bin " + lStream.str();
        _profhistPsi_pT[k]  = bookProfile1D(6+i+1, 1, j+1, hist_title[k][1]);
      }
    }
    
    string hist_title_Psi = "Psi(0.3 over R)";
    _profhistPsi = bookProfile1D(13, 1, 1, hist_title_Psi);
  }
  
  
  // Do the analysis
  void CDF_2005_S6217184::analyze(const Event& event) {
    Log log = getLog();
    log << Log::DEBUG << "Starting analyzing" << endl;
    
    // Analyse and print some info  
    const FastJets& jetpro = event.applyProjection(_jetsproj);
    
    log << Log::DEBUG << "Jet multiplicity before any pT cut = " << jetpro.getNJets() << endl;
    
    
    // Find vertex and check  that its z-component is < 60 cm from the nominal IP
    //const PVertex& pv = event.applyProjection(_vertexproj);
    /// @todo SEGV: either the HepMC event record is not filled properly or the F77-Wrapper functions are faulty
    /// @todo z- value assumed to be in mm, PYTHIA convention: dangerous!
    //if (fabs(pv.getPrimaryVertex().position().z()) < 600.0) {
    
    
    /// @todo Don't expose FastJet objects in Rivet analyses: the FastJets projection
    /// should convert them to Rivet 4-momentum classes (or similar).
    typedef vector<fastjet::PseudoJet> Jets;
    const Jets& jets = jetpro.getJetsPt();
    
    log << Log::DEBUG << "jetlist size = " << jets.size() << endl;
    
    int Njet = 0;
    bool jetcutpass = false;
    
    for (Jets::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {
      log << Log::DEBUG << "List item pT = " << jt->perp() << " E=" << jt->E() << " pz=" << jt->pz() << endl;
      if (jt->perp() > 37. && fabs(jt->rapidity())>0.1 && fabs(jt->rapidity())<0.7) jetcutpass = true;
      ++Njet;
      log << Log::DEBUG << "Jet pT =" << jt->perp() << " y=" << jt->rapidity() << " phi=" << jt->phi() << endl; 
    }
    
    if (jetcutpass) {
      
      const TotalVisibleMomentum& caloMissEt = event.applyProjection(_calmetproj);
      log << Log::DEBUG << "CaloMissEt.getMomentum().pT() = " << caloMissEt.getMomentum().pT() << endl;
      
      if (caloMissEt.getMomentum().pT()/sqrt(caloMissEt.getSET()) < 3.5) {
        
        FourMomentum jetaxis;
        _jetaxes.clear();
        for (Jets::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {
          // Only Central Calorimeter jets
          if (fabs(jt->rapidity()) < 1.1) {
            jetaxis.px(jt->px());
            jetaxis.py(jt->py());
            jetaxis.pz(jt->pz());
            jetaxis.E(jt->E());
            _jetaxes.push_back(jetaxis);
          }
        }
        if (_jetaxes.size()>0) { //determine jet shapes 
          const JetShape& jetShape = event.applyProjection(_jetshapeproj);
	  
          // Fill histograms.
          for (unsigned int jind=0; jind<_jetaxes.size(); ++jind) {
            for (int ipT=0; ipT<18; ++ipT) {
              if (_jetaxes[jind].pT() > _pTbins[ipT] && _jetaxes[jind].pT() <= _pTbins[ipT+1]) {
                _ShapeWeights[ipT] += event.weight(); 
                for (int rbin=0; rbin<jetShape.getNbins(); ++rbin) {
                  double rad_Rho = jetShape.getRmin() +(rbin+0.5)*jetShape.getInterval();
                  _profhistRho_pT[ipT]->fill(rad_Rho/_Rjet, jetShape.getDiffJetShape(jind, rbin), event.weight() );
                  double rad_Psi = jetShape.getRmin() +(rbin+1.0)*jetShape.getInterval();
                  _profhistPsi_pT[ipT]->fill(rad_Psi/_Rjet, jetShape.getIntJetShape(jind, rbin), event.weight() );
                }
                _profhistPsi->fill((_pTbins[ipT]+_pTbins[ipT+1])/2., jetShape.getPsi(jind), event.weight());
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
  void CDF_2005_S6217184::finalize() { 
    
  }
  
  
}
