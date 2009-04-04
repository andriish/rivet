// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2006_S6653332.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"


namespace Rivet {


 void CDF_2006_S6653332::init() {
   // Book histograms
   _sigmaBJet      = bookHistogram1D(1, 1, 1, "\\sigma(Z + b jet)");
   _ratioBJetToZ   = bookHistogram1D(2, 1, 1, "\\sigma(Z + b jet) / \\sigma(Z)");
   _ratioBJetToJet = bookHistogram1D(3, 1, 1, "\\sigma(Z + b jet) / \\sigma(Z + jet)");

  }  


  
  // Do the analysis
  void CDF_2006_S6653332::analyze(const Event& event) {
    // Check we have an l+l- pair that passes the kinematic cuts
     // Get the Z decay products (mu+mu- or e+e- pair)
    const InvMassFinalState& invMassFinalState = applyProjection<InvMassFinalState>(event, "INVFS");
    const ParticleVector&  ZDecayProducts =  invMassFinalState.particles();

    // make sure we have 2 Z decay products (mumu or ee) 
    if (ZDecayProducts.size() < 2) vetoEvent(event);

    _sumWeightsWithZ += event.weight();

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
	}
      }
    } // end loop around jets

    if(numJet > 0)    _sumWeightsWithZJet += event.weight();
    if(numBJet > 0) {
      _sigmaBJet->fill(1960.0,event.weight());
      _ratioBJetToZ->fill(1960.0,event.weight());
      _ratioBJetToJet->fill(1960.0,event.weight());

    }

  }    
    
  
  // Finalize
  void CDF_2006_S6653332::finalize() { 
    getLog() << Log::DEBUG << "Total sum of weights = " << sumOfWeights() << endl;
    getLog() << Log::DEBUG << "Sum of weights for Z production in mass range = " << _sumWeightsWithZ << endl;
    getLog() << Log::DEBUG << "Sum of weights for Z+jet production in mass range = " << _sumWeightsWithZJet << endl;

    _sigmaBJet->scale(crossSection()/sumOfWeights());
    _ratioBJetToZ->scale(1.0/_sumWeightsWithZ);
    _ratioBJetToJet->scale(1.0/_sumWeightsWithZJet);

  }
  
  
}
