// -*- C++ -*-
#include "Rivet/Analyses/D0_2007_S7075677.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  D0_2007_S7075677::D0_2007_S7075677()
  {
    // Run II Z rapidity
    setBeams(PROTON, ANTIPROTON);
    
    // All-inclusive final state
    FinalState fs;
    addProjection(fs, "FS");

    // Leading electrons
    LeadingParticlesFinalState lpfs(fs);
    lpfs.addParticleId(ELECTRON).addParticleId(POSITRON);
    addProjection(lpfs, "LeadingElectrons");

    IdentifiedFinalState ifs(fs);
    ifs.acceptId(PHOTON);
    addProjection(ifs, "PhotonFS");
  } 



  // Book histograms
  void D0_2007_S7075677::init() {
    _h_yZ = bookHistogram1D(1, 1, 1, "$1/\\sigma \\mathrm{d}\\sigma/\\mathrm{d}y(Z)$");
  }



  // Do the analysis 
  void D0_2007_S7075677::analyze(const Event & event) {
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
    // const ParticleVector allparts = applyProjection<FinalState>(event, "FS").particles();
    // const double HALO_RADIUS = 0.2;
    // foreach (const Particle& p, allparts) {
    //   if (p.pdgId() == PHOTON) {
    //     const double pho_eta = p.momentum().pseudorapidity();
    //     const double pho_phi = p.momentum().azimuthalAngle();
    //     /// NB. Not bothered about double-counting photons lying between the electrons
    //     if (deltaR(e1.pseudorapidity(), e1.azimuthalAngle(), pho_eta, pho_phi) < HALO_RADIUS) {
    //       e1 += p.momentum();
    //       continue;
    //     }
    //     if (deltaR(e2.pseudorapidity(), e2.azimuthalAngle(), pho_eta, pho_phi) < HALO_RADIUS) {
    //       e2 += p.momentum(); 
    //     }
    //   }
    // }
    const ParticleVector photons = applyProjection<FinalState>(event, "PhotonFS").particles();
    const double HALO_RADIUS = 0.2;
    foreach (const Particle& p, photons) {
      const double pho_eta = p.momentum().pseudorapidity();
      const double pho_phi = p.momentum().azimuthalAngle();
      /// NB. Not bothered about double-counting photons lying between the electrons
      if (deltaR(e1.pseudorapidity(), e1.azimuthalAngle(), pho_eta, pho_phi) < HALO_RADIUS) {
        e1 += p.momentum();
        continue;
      }
      if (deltaR(e2.pseudorapidity(), e2.azimuthalAngle(), pho_eta, pho_phi) < HALO_RADIUS) {
        e2 += p.momentum(); 
      }
    }

    // Calculate the Z momentum and cut on mass window [71,111]
    const FourMomentum Zmom = e1 + e2;
    if (Zmom.mass()/GeV < 71 || Zmom.mass()/GeV > 111) {
      getLog() << Log::DEBUG << "Electrons fall outside Z mass window" << endl;
      vetoEvent(event);
    }

    // Fill Z rapidity spectrum
    _h_yZ->fill(Zmom.rapidity(), weight);
  }



  // Finalize
  void D0_2007_S7075677::finalize() {
    normalize(_h_yZ);
  }


}
