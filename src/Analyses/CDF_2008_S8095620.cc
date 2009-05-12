// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2008_S8095620.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"

namespace Rivet {


  CDF_2008_S8095620::CDF_2008_S8095620()
    : _Rjet(0.7), _JetPtCut(20.), _JetEtaCut(1.5),  _sumWeightSelected(0.0)
  { 
    setBeams(PROTON, ANTIPROTON);
    
    const FinalState fs(-3.6, 3.6);
    addProjection(fs, "FS");
    
    // Create a final state with any e+e- or mu+mu- pair with 
    // invariant mass 76 -> 106 GeV and ET > 20 (Z decay products)
    std::vector<std::pair<long,long> > vids;
    vids.push_back(make_pair(ELECTRON, POSITRON));
    vids.push_back(make_pair(MUON, ANTIMUON));
    FinalState fs2(-3.6, 3.6);
    InvMassFinalState invfs(fs2, vids, 76*GeV, 106*GeV);
    addProjection(invfs, "INVFS");
    // Make a final state without the Z decay products for jet clustering
    VetoedFinalState vfs(fs);
    vfs.addVetoOnThisFinalState(invfs);
    addProjection(vfs, "VFS");
    addProjection(FastJets(vfs, FastJets::CDFMIDPOINT, 0.7), "Jets");
  }


  void CDF_2008_S8095620::init() {
     // Book histograms
    _dSdET    = bookHistogram1D(1, 1, 1, "1/\\sigma(Z) \\times d\\sigma(Z + b jet)/dE_{T}^{jet}",
				"$E_{T}^{jet}$ [GeV]","1/\\sigma(Z) \\times d\\sigma(Z + b jet)/dE_{T}^{jet}");
    _dSdETA   = bookHistogram1D(2, 1, 1, "1/\\sigma(Z) \\times d\\sigma(Z + b jet)/d\\eta^{jet}",
				"$\\eta^{jet}$","1/\\sigma(Z) \\times d\\sigma(Z + b jet)/\\eta^{jet}");
    _dSdNJet  = bookHistogram1D(3, 1, 1, "1/\\sigma(Z) \\times d\\sigma(Z + b jet)/dN^{jet}",
				"Number of jets","1/\\sigma(Z) \\times d\\sigma(Z + b jet)/dN^{jet}");
    _dSdNbJet = bookHistogram1D(4, 1, 1, "1/\\sigma(Z) \\times d\\sigma(Z + b jet)/dN^{b jet} ",
				"Number of b jets","1/\\sigma(Z) \\times d\\sigma(Z + b jet)/dN^{b jet}");
    _dSdZpT   = bookHistogram1D(5, 1, 1, "1/\\sigma(Z) \\times d\\sigma(Z + b jet)/dp_T{}^{Z}",
				"Z $p_{T}$ [GeV]","1/\\sigma(Z) \\times d\\sigma(Z + b jet)/dp_{T}^{Z}");
  }  


  
  // Do the analysis
  void CDF_2008_S8095620::analyze(const Event& event) {
    // Check we have an l+l- pair that passes the kinematic cuts
     // Get the Z decay products (mu+mu- or e+e- pair)
    const InvMassFinalState& invMassFinalState = applyProjection<InvMassFinalState>(event, "INVFS");
    const ParticleVector&  ZDecayProducts =  invMassFinalState.particles();

    // make sure we have 2 Z decay products (mumu or ee) 
    if (ZDecayProducts.size() < 2) vetoEvent(event);
    _sumWeightSelected += event.weight(); 
    // @todo: write out a warning if there are more than two decay products
    FourMomentum Zmom = ZDecayProducts[0].momentum() +  ZDecayProducts[1].momentum();

    // Put all b-quarks in a vector
    ParticleVector bquarks;
    /// @todo Provide nicer looping
    for (GenEvent::particle_const_iterator p = event.genEvent().particles_begin(); 
         p != event.genEvent().particles_end(); ++p) {
      if ( fabs((*p)->pdg_id()) == BQUARK ) {
        bquarks.push_back(Particle(**p));
      }
    }
    
    // Get jets 
    const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
    getLog() << Log::DEBUG << "Jet multiplicity before any pT cut = " << jetpro.size() << endl;
    
    const PseudoJets& jets = jetpro.pseudoJetsByPt();
    getLog() << Log::DEBUG << "jetlist size = " << jets.size() << endl;

    int numBJet = 0;
    int numJet  = 0;
    // for each b-jet plot the ET and the eta of the jet, normalise to the total cross section at the end
    // for each event plot N jet and pT(Z), normalise to the total cross section at the end 
    for (PseudoJets::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {
      // select jets that pass the kinematic cuts
      if (jt->perp() > _JetPtCut && fabs(jt->rapidity()) <= _JetEtaCut) {
	numJet++;
	// does the jet contain a b-quark?
	bool bjet = false;
	foreach (const Particle& bquark,  bquarks) {
	  if (deltaR(jt->rapidity(), jt->phi(), bquark.momentum().rapidity(),bquark.momentum().azimuthalAngle()) <= _Rjet) {
	    bjet = true;
	    break;
	  }
	} // end loop around b-jets
	if (bjet) {
          numBJet++;
	  _dSdET->fill(jt->perp(),event.weight()); 
	  _dSdETA->fill(jt->rapidity(),event.weight()); 
	}
      }
    } // end loop around jets

    if(numJet > 0) _dSdNJet->fill(numJet,event.weight());
    if(numBJet > 0) {
      _dSdNbJet->fill(numBJet,event.weight());
      _dSdZpT->fill(Zmom.pT(),event.weight()); 
    }    
  }

  

  // Finalize
  void CDF_2008_S8095620::finalize() {  
    // normalise histograms
    // scale by 1 / the sum-of-weights of events that pass the Z cuts
    // since the cross sections are normalized to the inclusive
    // Z cross sections. 
    double Scale = 1.0;
    if (_sumWeightSelected != 0.0) Scale = 1.0/_sumWeightSelected;
    _dSdET->scale(Scale);
    _dSdETA->scale(Scale);
    _dSdNJet->scale(Scale);
    _dSdNbJet->scale(Scale);
    _dSdZpT->scale(Scale);
  }
}
