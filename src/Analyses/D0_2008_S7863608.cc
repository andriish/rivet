// -*- C++ -*-
#include "Rivet/Analyses/D0_2008_S7863608.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  D0_2008_S7863608::D0_2008_S7863608()
  {
    setBeams(PROTON, ANTIPROTON);
    
    /// @todo Use cross-section from generator
    //setNeedsCrossSection(true);

    //full final state
    FinalState fs(-5.0, 5.0);
    addProjection(fs, "FS");
 
    //leading muons in tracking acceptance
    LeadingParticlesFinalState lpfs(fs, -1.7, 1.7, 15.0*GeV);
    lpfs.addParticleId(MUON).addParticleId(ANTIMUON);
    addProjection(lpfs, "LeadingMuons");

    // Vetoed fs for jets
    VetoedFinalState vfs(fs);
    // Veto the muons from Z decay  
    vfs.addVetoOnThisFinalState(lpfs);
    addProjection(vfs, "VFS");
  } 



  // Book histograms
  void D0_2008_S7863608::init() {

    /// @todo Dividing through by measured Z cross-section would be nice...
    _h_jet_pT_cross_section = bookHistogram1D(1, 1, 1, "Differential cross section in leading jet $p_\\perp$");
    _h_jet_y_cross_section = bookHistogram1D(2, 1, 1, "Differential cross section in leading jet rapidity");
    _h_Z_pT_cross_section = bookHistogram1D(3, 1, 1, "Differential cross section in Z/$\\gamma*$ $p_\\perp$");
    _h_Z_y_cross_section = bookHistogram1D(4, 1, 1, "Differential cross section in Z/$\\gamma*$ rapidity");
    _h_total_cross_section = bookHistogram1D(5, 1, 1, "Total Z + jet cross section");
    
  }



  // Do the analysis 
  void D0_2008_S7863608::analyze(const Event & event) {
    double weight = event.weight();

    // Skip if the event is empty
    const FinalState& fs = applyProjection<FinalState>(event, "FS");
    if (fs.isEmpty()) {
      getLog() << Log::DEBUG << "Skipping event " << event.genEvent().event_number()
               << " because no final state pair found " << endl;
      vetoEvent(event);
    }


    // Find the Z candidates
    const FinalState & muonfs = applyProjection<FinalState>(event, "LeadingMuons");
    // If there are no muons in the FinalState, skip the event
    if (muonfs.particles().size() != 2) {
      getLog() << Log::DEBUG << "Skipping event " << event.genEvent().event_number()
               << " because no muon pair found " << endl;
      vetoEvent(event);
    }
    
    // Now build the jets on a FS without the muons from the Z
    const ParticleVector vparts = applyProjection<FinalState>(event, "VFS").particles();
    ParticleVector jetparts;
    foreach (const Particle& p, vparts) {
      if (p.pdgId() == PHOTON) {
        bool copy = true;
        foreach (const Particle& mu, muonfs.particles()) {
          if (deltaR(mu.momentum().pseudorapidity(), mu.momentum().azimuthalAngle(),
                     p.momentum().pseudorapidity(), p.momentum().azimuthalAngle()) < 0.2) {
            copy = false;
            break;
          }
        }
        if (!copy) {
          getLog() << Log::DEBUG << "Excluding photon from muon" << endl;
          continue;
        }
      }
      jetparts.push_back(p);
    }

    /// @todo Allow proj creation w/o FS as ctor arg, so that calc can be used more easily.
    D0ILConeJets jetpro(muonfs, 0.5); //< @todo The 'muonfs' arg makes no sense!
    jetpro.calc(jetparts);

    // Take the leading jet with pt > 20, |y| < 2.8:
    /// @todo Make this neater, using the JetAlg interface and the built-in sorting
    const Jets& jets = jetpro.jets();
    Jets jets_cut;
    foreach (const Jet& j, jets) {
      if (j.momentum().pT()/GeV > 20 && fabs(j.momentum().rapidity()) < 2.8) {
        jets_cut.push_back(j);
      }
    }
    getLog() << Log::DEBUG << "Num jets above 20 GeV = " << jets_cut.size() << endl;

    // Return if there are no jets:
    if(jets_cut.size()<1) {
      getLog() << Log::DEBUG << "Skipping event " << event.genEvent().event_number()
               << " because no jets pass cuts " << endl;
      vetoEvent(event);
    }

    // Sort by pT:
    sort(jets_cut.begin(), jets_cut.end(), cmpJetsByPt);

    // Calculate the Z pT, rapidity:
    const ParticleVector muons = muonfs.particles();
    const FourMomentum Zmom = muons[0].momentum() + muons[1].momentum();

    if (Zmom.mass()/GeV < 65 || Zmom.mass()/GeV > 115) {
      getLog() << Log::DEBUG << "Skipping event " << event.genEvent().event_number()
               << " because failing Z mass window " << endl;
      vetoEvent(event);
    }

    // cut on Delta R between jet and muons
    foreach (const Jet& j, jets_cut) {
      foreach (const Particle& mu, muons) {
        if (deltaR(mu.momentum().pseudorapidity(), mu.momentum().azimuthalAngle(),
                   j.momentum().pseudorapidity(), j.momentum().azimuthalAngle()) < 0.5) {
          vetoEvent(event);
        }
      }
    }

    // In jet pT
    _h_jet_pT_cross_section->fill( jets_cut[0].momentum().pT(), weight);
    _h_jet_y_cross_section->fill( fabs(jets_cut[0].momentum().rapidity()), weight);

    // In Z pT
    _h_Z_pT_cross_section->fill(Zmom.pT(), weight);
    _h_Z_y_cross_section->fill(fabs(Zmom.rapidity()), weight);

    // _h_total_cross_section = bookHistogram1D     _crossSectionRatio->fill(1, weight);

  }



  // Finalize
  void D0_2008_S7863608::finalize() {
    /// @todo Use the generator cross-section
    //_h_total_cross_section->fill(crossSection());
    normalize(_h_jet_pT_cross_section, 18.7);
    normalize(_h_jet_y_cross_section, 18.7);
    normalize(_h_Z_pT_cross_section, 18.7);
    normalize(_h_Z_y_cross_section, 18.7);
  }

}
