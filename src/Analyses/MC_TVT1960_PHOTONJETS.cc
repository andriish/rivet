// -*- C++ -*-
#include "Rivet/Analyses/MC_TVT1960_PHOTONJETS.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  MC_TVT1960_PHOTONJETS::MC_TVT1960_PHOTONJETS()
    : MC_JetAnalysis("MC_TVT1960_PHOTONJETS", 1960.0, 4, "Jets")
  {
    setBeams(PROTON, ANTIPROTON);
    setNeedsCrossSection(true);
    
    // General FS
    FinalState fs(-5.0, 5.0);
    addProjection(fs, "FS");
 
    // Get leading photon
    LeadingParticlesFinalState photonfs(fs, -1.0, 1.0);
    photonfs.addParticleId(PHOTON);
    addProjection(photonfs, "LeadingPhoton");

    // FS for jets excludes the leading photon
    VetoedFinalState vfs(fs);
    vfs.addVetoOnThisFinalState(photonfs);
    addProjection(vfs, "JetFS");
    FastJets jetpro(vfs, FastJets::KT, 0.7, 20.0*GeV);
    addProjection(jetpro, "Jets");
  }



  // Book histograms
  void MC_TVT1960_PHOTONJETS::init() {
    _h_photon_pT = bookHistogram1D("photon_pT", 100, 0.0, 500.0);
    _h_photon_y = bookHistogram1D("photon_y", 40, -4.0, 4.0);
    _h_photon_jet1_deta = bookHistogram1D("photon_jet1_deta", 50, -5.0, 5.0);
    _h_photon_jet1_dR = bookHistogram1D("photon_jet1_dR", 25, 0.5, 7.0);
    
    MC_JetAnalysis::init();
  }



  // Do the analysis 
  void MC_TVT1960_PHOTONJETS::analyze(const Event & e) {
    double weight = e.weight();
    
    // Get the photon
    const ParticleVector photons = applyProjection<FinalState>(e, "LeadingPhoton").particles();
    if (photons.size() != 1) {
      vetoEvent;
    }
    const FourMomentum photon = photons.front().momentum();
    if (photon.pT()/GeV < 30) {
      getLog() << Log::DEBUG << "Leading photon has pT < 30 GeV: " << photon.pT()/GeV << endl;
      vetoEvent;
    }

    // Get all charged particles
    const FinalState& fs = applyProjection<FinalState>(e, "JetFS");
    if (fs.isEmpty()) {
      vetoEvent;
    }

    // Isolate photon by ensuring that a 0.4 cone around it contains less than 7% of the photon's energy
    const double egamma = photon.E();
    double econe = 0.0;
    foreach (const Particle& p, fs.particles()) {
      if (deltaR(photon, p.momentum()) < 0.4) {
        econe += p.momentum().E();
        // Veto as soon as E_cone gets larger
        if (econe/egamma > 0.07) {
          vetoEvent;
        }
      }
    }
    
    _h_photon_pT->fill(photon.pT(),weight);
    _h_photon_y->fill(photon.rapidity(),weight);
    
    const FastJets& jetpro = applyProjection<FastJets>(e, "Jets");
    const Jets& jets = jetpro.jetsByPt(20.0*GeV);
    if (jets.size()>0) {
      _h_photon_jet1_deta->fill(photon.eta()-jets[0].momentum().eta(), weight);
      _h_photon_jet1_dR->fill(deltaR(photon, jets[0].momentum()), weight);
    }

    MC_JetAnalysis::analyze(e);
  }


  // Finalize
  void MC_TVT1960_PHOTONJETS::finalize() {
    scale(_h_photon_pT, crossSection()/sumOfWeights());
    scale(_h_photon_y, crossSection()/sumOfWeights());
    scale(_h_photon_jet1_deta, crossSection()/sumOfWeights());
    scale(_h_photon_jet1_dR, crossSection()/sumOfWeights());
    
    MC_JetAnalysis::finalize();
  }

}
