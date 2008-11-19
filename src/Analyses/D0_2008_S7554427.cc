// -*- C++ -*-
#include "Rivet/Analyses/D0_2008_S7554427.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  D0_2008_S7554427::D0_2008_S7554427()
  {
    setBeams(PROTON, ANTIPROTON);
    
    /// @todo Use cross-section from generator
    //setNeedsCrossSection(true);

    // All-inclusive final state
    FinalState fs;
    addProjection(fs, "FS");

    // Leading electrons
    LeadingParticlesFinalState lpfs(fs);
    lpfs.addParticleId(ELECTRON).addParticleId(POSITRON);
    addProjection(lpfs, "LeadingElectrons");
  } 



  // Book histograms
  void D0_2008_S7554427::init() {
    _h_ZpT         = bookHistogram1D(1, 1, 1, "$1/\\sigma \\mathrm{d}\\sigma/\\mathrm{d}p_\\perp(Z)$");
    _h_forward_ZpT = bookHistogram1D(3, 1, 1, "$1/\\sigma \\mathrm{d}\\sigma/\\mathrm{d}p_\\perp(Z)$ (forward region only)");
  }



  // Do the analysis 
  void D0_2008_S7554427::analyze(const Event & event) {
    const double weight = event.weight();

    // Find the Z candidates
    const FinalState& efs = applyProjection<FinalState>(event, "LeadingElectrons");
    // If there are no electrons in the FinalState, skip the event
    if (efs.particles().size() != 2) {
      getLog() << Log::DEBUG << "No electron pair found " << endl;
      vetoEvent(event);
    }
    FourMomentum e1 = efs.particles()[0].momentum();
    FourMomentum e2 = efs.particles()[1].momentum();
    
    // Add back in photons that could have radiated from the Z decay products
    const ParticleVector allparts = applyProjection<FinalState>(event, "FS").particles();
    foreach (const Particle& p, allparts) {
      if (p.pdgId() == PHOTON) {
        const double pho_eta = p.momentum().pseudorapidity();
        const double pho_phi = p.momentum().azimuthalAngle();
        /// @todo Need to be super-careful about photons lying between the electrons? (no)
        if (deltaR(e1.pseudorapidity(), e1.azimuthalAngle(), pho_eta, pho_phi) < 0.2) {
          e1 += p.momentum();
          continue;
        }
        if (deltaR(e2.pseudorapidity(), e2.azimuthalAngle(), pho_eta, pho_phi) < 0.2) {
          e2 += p.momentum(); 
        }
      }
    }

    // Calculate the Z momentum and cut on mass window [40,200]
    const FourMomentum Zmom = e1 + e2;
    if (Zmom.mass()/GeV < 40 || Zmom.mass()/GeV > 200) {
      getLog() << Log::DEBUG << "Electrons fall outside Z mass window" << endl;
      vetoEvent(event);
    }

    // Fill Z pT spectra
    _h_ZpT->fill(Zmom.pT(), weight);
    if (fabs(Zmom.rapidity()) > 2) {
      _h_forward_ZpT->fill(Zmom.pT(), weight);
    }
  }



  // Finalize
  void D0_2008_S7554427::finalize() {
    normalize(_h_ZpT);
    normalize(_h_forward_ZpT);
  }

}
