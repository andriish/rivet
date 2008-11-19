// -*- C++ -*-
#include "Rivet/Analyses/D0_2008_S7719523.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  D0_2008_S7719523::D0_2008_S7719523() : _events(0) 
  {
    setBeams(PROTON, ANTIPROTON);
    
    /// @todo Use cross-section from generator
    //setNeedsCrossSection(true);

    // Forward and central particles for 
    FinalState fs(-5.0, 5.0);
    addProjection(fs, "FS");
 
    // Get leading photon
    /// @todo Repeating the cuts is a bit silly...

    LeadingParticlesFinalState photonfs(fs, -1.0, 1.0);
    photonfs.addParticleId(PHOTON);
    addProjection(photonfs, "LeadingPhoton");

    /// @todo Remove leading photon from jet tracks

    // //Vetoed fs for jets
    // VetoedFinalState vfs(fs);
    // // veto the muons from Z decay  
    // vfs.addVetoOnThisFinalState(lpfs);
    // addProjection(vfs, "VFS");

    /// @todo Remove extra photons or incorporate into jets?
    /// @todo Veto neutrinos
  } 



  // Book histograms
  void D0_2008_S7719523::init() {
    _h_central_same_cross_section = 
      bookHistogram1D(1, 1, 1, "$\\mathrm{d}\\sigma/\\mathrm{d}p_\\perp(\\gamma_\\text{lead})$ (central jets, same-sign rapidity)");
    _h_central_opp_cross_section  = 
      bookHistogram1D(2, 1, 1, "$\\mathrm{d}\\sigma/\\mathrm{d}p_\\perp(\\gamma_\\text{lead})$ (central jets, opp-sign rapidity)");
    _h_forward_same_cross_section = 
      bookHistogram1D(3, 1, 1, "$\\mathrm{d}\\sigma/\\mathrm{d}p_\\perp(\\gamma_\\text{lead})$ (forward jets, same-sign rapidity)");
    _h_forward_opp_cross_section  = 
      bookHistogram1D(4, 1, 1, "$\\mathrm{d}\\sigma/\\mathrm{d}p_\\perp(\\gamma_\\text{lead})$ (forward jets, opp-sign rapidity)"); 
  }



  // Do the analysis 
  void D0_2008_S7719523::analyze(const Event& event) {
    const double weight = event.weight();

    // Get the photon
    const FinalState& photonfs = applyProjection<FinalState>(event, "LeadingPhoton");
    if (photonfs.particles().size() != 1) {
      getLog() << Log::DEBUG << "No photon found" << endl;
      vetoEvent(event);
    }
    const FourMomentum photon = photonfs.particles().front().momentum();
    if (photon.pT()/GeV < 30) {
      getLog() << Log::DEBUG << "Leading photon has pT < 30 GeV: " << photon.pT()/GeV << endl;
      vetoEvent(event);
    }

    /// @todo Isolate photon by counting everything in a 0.4 cone around it and
    ///       requiring that it has less than ???% of the photon's energy

    // Get all charged particles
    const FinalState& fs = applyProjection<FinalState>(event, "FS");
    if (fs.isEmpty()) {
      vetoEvent(event);
    }
    /// @todo Allow proj creation w/o FS as ctor arg, so that calc can be used more easily.
    D0ILConeJets jetpro(fs); //< @todo This arg makes no sense!
    jetpro.calc(fs.particles());
    Jets isolated_jets;
    foreach (const Jet& j, jetpro.jets()) {
      const FourMomentum pjet = j.momentum();
      const double dr = deltaR(photon.pseudorapidity(), photon.azimuthalAngle(),
                               pjet.pseudorapidity(), pjet.azimuthalAngle());
      if (dr > 0.7 && pjet.pT()/GeV > 15) {
        isolated_jets.push_back(j);
      }
    }
    
    getLog() << Log::DEBUG << "Num jets after isolation and pT cuts = " << isolated_jets.size() << endl;
    if (isolated_jets.empty()) {
      getLog() << Log::DEBUG << "No jets pass cuts" << endl;
      vetoEvent(event);
    }

    // Sort by pT and get leading jet
    sort(isolated_jets.begin(), isolated_jets.end(), cmpJetsByPt);
    const FourMomentum leadingJet = isolated_jets.front().momentum();
    int photon_jet_sign = sign( leadingJet.rapidity() * photon.rapidity() );

    /// @todo What if leading jet is in [0.8, 1.5]?

    if (fabs(leadingJet.rapidity()) < 0.8) {
      if (photon_jet_sign >= 1) {
        _h_central_same_cross_section->fill(photon.pT(), weight);
      } else {
        _h_central_opp_cross_section->fill(photon.pT(), weight);
      }
    } else if (inRange( fabs(leadingJet.rapidity()), 1.5, 2.5)) {
      if (photon_jet_sign >= 1) {
        _h_forward_same_cross_section->fill(photon.pT(), weight);
      } else {
        _h_forward_opp_cross_section->fill(photon.pT(), weight); 
      }
    }

    /// @todo Cross-section ratios (6 plots)

  }



  // Finalize
  void D0_2008_S7719523::finalize() {
    /// @todo Use the generator cross-section
    //_h_total_cross_section->fill(crossSection());
    normalize(_h_central_same_cross_section, 347.4);
    normalize(_h_central_opp_cross_section,  281.8);
    normalize(_h_forward_same_cross_section, 164.8);
    normalize(_h_forward_opp_cross_section,   81.5);
  }

}
