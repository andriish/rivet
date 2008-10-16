// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2008_S7782535.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
// @todo : test with Pythia
namespace Rivet {

  void CDF_2008_S7782535::init() {
    _pTbins.push_back(52.);
    _pTbins.push_back(80.);
    _pTbins.push_back(104.);
    _pTbins.push_back(142.);
    _pTbins.push_back(300.);
     // Book histograms
    for (int i = 0; i < _NpTbins ; i++) {
       stringstream title;
       title << "Integral jet shape Psi," << _pTbins[i] << " < pT < "<< _pTbins[i+1]; 
         _Psi_pT[i] = bookProfile1D(i+1,2,1,title.str());
    }
      _OneMinusPsi_vs_pT = bookDataPointSet(5,1,1,"1 - Psi vs Jet pT");
  }  

  
  // Do the analysis
  void CDF_2008_S7782535::analyze(const Event& event) {
    Log log = getLog();
    log << Log::DEBUG << "Starting analyzing" << endl;

    // put all b-quarks in a vector
    ParticleVector bquarks;
    for (GenEvent::particle_const_iterator p = event.genEvent().particles_begin();
	 p != event.genEvent().particles_end(); ++p) 
      if ( fabs((*p)->pdg_id())  == 5 ) bquarks.push_back(Particle(**p));

    if (!bquarks.size()) { 
      log << Log::DEBUG << "No b-quarks, exiting" << endl;
      vetoEvent(event);
    }
    // Get final state particles in event  
    ////    const FinalState& part = applyProjection<FinalState>(event, "FS");
    ////    const ParticleVector& particles =  part.particles();

    // Get jets 
    const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
    log << Log::DEBUG << "Jet multiplicity before any pT cut = " << jetpro.getNumJets() << endl;
    
    /// @todo Don't expose FastJet objects in Rivet analyses: the FastJets projection
    /// should convert them to Rivet 4-momentum classes (or similar).
    const PseudoJets& jets = jetpro.pseudoJetsByPt();
    log << Log::DEBUG << "jetlist size = " << jets.size() << endl;
    // Determine the central jet axes
    FourMomentum jetaxis;
    _jetaxes.clear();
    for (PseudoJets::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {
      // Only Central Calorimeter jets
      /// @todo Declare this cut
      //      cout<<jt->perp() << " " << jt->rapidity() << endl;
      if (jt->perp() > _pTbins[0] && fabs(jt->rapidity()) <= 0.7) {
        jetaxis.px(jt->px());
        jetaxis.py(jt->py());
        jetaxis.pz(jt->pz());
        jetaxis.E(jt->E());
        _jetaxes.push_back(jetaxis);
      }
    }
    // Determine jet shapes
    if (_jetaxes.empty())  { 
      log << Log::DEBUG << "No jet axes" << endl;
      vetoEvent(event);
    }
      
    const JetShape& jetShape = applyProjection<JetShape>(event, "JetShape");
    for (size_t jind = 0; jind < _jetaxes.size(); ++jind) { 
      bool bjet = false;
      for (ParticleVector::const_iterator bquark =  bquarks.begin();
	   bquark != bquarks.end() && !bjet; ++bquark) {
	if (deltaR( _jetaxes[jind].rapidity(), _jetaxes[jind].azimuthalAngle(), bquark->momentum().rapidity(), bquark->momentum().azimuthalAngle()) <= _Rjet ) bjet=true;
      } 
      if(bjet) {	
	// put jet in correct pT bin
	int jet_pt_bin = -1;
	if      (_jetaxes[jind].pT() > _pTbins[0] && _jetaxes[jind].pT() <= _pTbins[1]) jet_pt_bin = 0;
	else if (_jetaxes[jind].pT() > _pTbins[1] && _jetaxes[jind].pT() <= _pTbins[2]) jet_pt_bin = 1;
	else if (_jetaxes[jind].pT() > _pTbins[2] && _jetaxes[jind].pT() <= _pTbins[3]) jet_pt_bin = 2;
	else if (_jetaxes[jind].pT() > _pTbins[3] && _jetaxes[jind].pT() <= _pTbins[4]) jet_pt_bin = 3;
	if (jet_pt_bin > -1) {
	  // fill each entry in profile
	  for (size_t rbin = 0; rbin < jetShape.getNbins(); ++rbin) {
	    const double rad_Psi = jetShape.getRmin() +(rbin+1.0)*jetShape.getInterval();
	    _Psi_pT[jet_pt_bin]->fill(rad_Psi/_Rjet, jetShape.getIntJetShape(jind, rbin), event.weight() );
	  }
	} // end valid jet_pt_bin
      } // end bjet
    } // end loop round jets
  
}

  

  // Finalize
  void CDF_2008_S7782535::finalize() {  
    std::vector<double> x, y, ex, ey;
    for (unsigned int i = 0; i < _pTbins.size()-1; i++) {
      x.push_back((_pTbins[i]+_pTbins[i+1])/2.);
      ex.push_back(0.);
      // get entry for  rad_Psi = 0.2 bin
      y.push_back(1.0 - _Psi_pT[i]->binHeight(1));
      ey.push_back(_Psi_pT[i]->binError(1)); 
      
    }
    _OneMinusPsi_vs_pT->setCoordinate(0,x,ex);
    _OneMinusPsi_vs_pT->setCoordinate(1,y,ey);

  }
}
