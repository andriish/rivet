// -*- C++ -*-
// CDF Z pT analysis

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2000_S4155203.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/LossyFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  CDF_2000_S4155203::CDF_2000_S4155203() { 
    setBeams(PROTON, ANTIPROTON);
    const FinalState fs(-4.2, 4.2, 15*GeV);

    // Get photons  
    IdentifiedFinalState ifs(fs);
    ifs.acceptId(PHOTON);
    addProjection(ifs, "PhotonFS");

    // Leading electrons
    // LeadingParticlesFinalState lpfs(fs);
    // lpfs.addParticleId(ELECTRON).addParticleId(POSITRON);
    // addProjection(lpfs, "LeadingElectrons");  
    addProjection(ChargedLeptons(fs), "CL");
  }



  // Book histograms
  void CDF_2000_S4155203::init() {
    _hist_zpt = bookHistogram1D(1, 1, 1, "$p_\\perp$ of Z boson in $\\Pelectron \\Ppositron$ decays",
                                "$p_\\perp(\\PZ)$ / GeV",
                                "$\\d{\\sigma}/\\d{p_\\perp(\\PZ)}$"); //< @todo Cross-section units
  }


  // Do the analysis
  void CDF_2000_S4155203::analyze(const Event& e) {
    const double weight = e.weight();

    // Get the leptons
    const ParticleVector& leptons = applyProjection<ChargedLeptons>(e, "CL").chargedLeptons();
    ParticleVector electrons;
    electrons.clear();

    bool has_central_electron = false;
    foreach (const Particle& l, leptons) {
      // Only take electrons
      if (abs(l.pdgId()) != ELECTRON) continue;
      
      // Check for electrons in the fiducial regions
      if (fabs(l.momentum().pseudorapidity()) < 1.1 &&
          l.momentum().pT()/GeV > 20) {
        has_central_electron = true;
        electrons.push_back(l);
      }
      else if (fabs(l.momentum().pseudorapidity()) > 1.1 &&
               fabs(l.momentum().pseudorapidity()) < 4.2 &&
               l.momentum().pT()/GeV > 15) {
        electrons.push_back(l);
      }
    }

    // One electron needs to be in the central region
    if (! has_central_electron) {
      vetoEvent(e);
    }

    // We want exactly two electrons with opposite charge
    getLog() << Log::DEBUG << "Electron multiplicity = " << electrons.size() << endl;
    if (electrons.size() != 2 || electrons[0].pdgId() != -electrons[1].pdgId() ) {
      vetoEvent(e);
    }

    
    // Add back in photons that could have radiated from the Z decay products
    FourMomentum e1 = electrons[0].momentum();
    FourMomentum e2 = electrons[1].momentum();
    const ParticleVector photons = applyProjection<FinalState>(e, "PhotonFS").particles();
    const double HALO_RADIUS = 0.2;
    foreach (const Particle& p, photons) {
      const double pho_eta = p.momentum().pseudorapidity();
      const double pho_phi = p.momentum().azimuthalAngle();
      // NB. Not bothered about double-counting photons lying between the electrons
      if (deltaR(e1.pseudorapidity(), e1.azimuthalAngle(), pho_eta, pho_phi) < HALO_RADIUS) {
        e1 += p.momentum();
        continue;
      }
      if (deltaR(e2.pseudorapidity(), e2.azimuthalAngle(), pho_eta, pho_phi) < HALO_RADIUS) {
        e2 += p.momentum(); 
      }
    }


    FourVector dilepton = electrons[0].momentum() + electrons[1].momentum();
    const double mZ = mass(dilepton);
    const double pTZ = pT(dilepton);

    // Lepton pair should have an invariant mass between 66 and 116
    if (mZ/GeV < 66.0 || mZ/GeV > 116.0) {
      vetoEvent(e);
    }

    getLog() << Log::DEBUG << "Dilepton mass = " << mZ/GeV << " GeV"  << endl;
    getLog() << Log::DEBUG << "Dilepton pT   = " << pTZ/GeV << " GeV" << endl;

    // Fill the histograms
    _hist_zpt->fill(pTZ/GeV, weight);
  }


  void CDF_2000_S4155203::finalize() {
    // Normalize to the experimental cross-section
    /// @todo Get norm from generator cross-section
    normalize(_hist_zpt, 247.4);
  }


}
