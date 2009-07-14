// -*- C++ -*-
#include "Rivet/Analyses/D0_2006_S6438750.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  D0_2006_S6438750::D0_2006_S6438750()
    : Analysis("D0_2006_S6438750")
  {
    setBeams(PROTON, ANTIPROTON);
    
    /// @todo Use cross-section from generator
    setNeedsCrossSection(true);

    // General FS for photon isolation
    FinalState fs(-1.5, 1.5);
    addProjection(fs, "AllFS");

    // Get leading photon
    LeadingParticlesFinalState photonfs(fs, -1.0, 1.0);
    photonfs.addParticleId(PHOTON);
    addProjection(photonfs, "LeadingPhoton");
  } 



  // Book histograms
  void D0_2006_S6438750::init() {
    _h_pTgamma = bookHistogram1D(1, 1, 1);
  }



  // Do the analysis 
  void D0_2006_S6438750::analyze(const Event& event) {

    // Get the photon
    const FinalState& photonfs = applyProjection<FinalState>(event, "LeadingPhoton");
    if (photonfs.particles().size() != 1) {
      getLog() << Log::DEBUG << "No photon found" << endl;
      vetoEvent;
    }
    const FourMomentum photon = photonfs.particles().front().momentum();
    if (photon.pT()/GeV < 23) {
      getLog() << Log::DEBUG << "Leading photon has pT < 23 GeV: " << photon.pT()/GeV << endl;
      vetoEvent;
    }

    // Get all other particles
    const FinalState& fs = applyProjection<FinalState>(event, "AllFS");
    if (fs.isEmpty()) {
      vetoEvent;
    }

    // Isolate photon by ensuring that a 0.4 cone around it contains less than 7% of the photon's energy
    const double egamma = photon.E();
    // Energy inside R = 0.2
    double econe_02 = 0.0;
    // Energy between R = [0.2, 0.4]
    double econe_02_04 = 0.0;
    foreach (const Particle& p, fs.particles()) {
      const double dr = deltaR(photon.pseudorapidity(), photon.azimuthalAngle(),
                               p.momentum().pseudorapidity(), p.momentum().azimuthalAngle());
      if (dr < 0.2) {
        econe_02 += p.momentum().E();
      } else if (dr < 0.4) {
        econe_02_04 += p.momentum().E();
      }
    }
    // Veto if outer hollow cone contains more than 10% of the energy of the inner cone
    // or if the non-photon energy in the inner cone exceeds 5% of the photon energy.
    if (econe_02_04/econe_02 > 0.1 || (econe_02-egamma)/egamma > 0.05) {
      getLog() << Log::DEBUG << "Vetoing event because photon is insufficiently isolated" << endl;
      vetoEvent;
    }

    // Veto if leading jet is outside plotted rapidity regions
    const double eta_gamma = fabs(photon.pseudorapidity());
    if (eta_gamma > 0.9) {
      getLog() << Log::DEBUG << "Leading photon falls outside acceptance range; "
               << "|eta_gamma| = " << eta_gamma << endl;
      vetoEvent;
    }
    
    // Fill histo
    const double weight = event.weight();
    _h_pTgamma->fill(photon.pT(), weight); 
  }



  // Finalize
  void D0_2006_S6438750::finalize() {
    /// @todo Generator cross-section from Pythia gives ~7500, vs. expected 2988!
    //normalize(_h_pTgamma, 2988.4869);

    const double lumi_gen = sumOfWeights()/crossSection();
    // Divide by effective lumi, plus rapidity bin width of 1.8
    scale(_h_pTgamma, 1/lumi_gen * 1/1.8);
  }

}
