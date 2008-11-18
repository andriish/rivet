// -*- C++ -*-
#include "Rivet/Analyses/D0_2008_S6879055.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {

  D0_2008_S6879055::D0_2008_S6879055() : _events(0) 
  {
    setBeams(PROTON, ANTIPROTON);
    //full final state
    FinalState fs(-5.0, 5.0);
    addProjection(fs, "FS");
    //leading electrons in tracking acceptance
    LeadingParticlesFinalState lpfs(fs, -1.1, 1.1, 25);
    lpfs.addParticleId(11).addParticleId(-11);
    //addProjection(lpfs, "LeadingElectronsFS");
    //InvMass selection
    InvMassFinalState electronsFromZ(lpfs, make_pair(11, -11), 75, 105, -1.1, 1.1, 25);
    addProjection(electronsFromZ,"ElectronsFromZ");
    //Vetoed fs for jets
    VetoedFinalState vfs(fs);
    // Add particle/antiparticle vetoing: 12=nu_e, 14=nu_mu, 16=nu_tau
    vfs.addVetoPairId(12).addVetoPairId(14).addVetoPairId(16);
    // veto the electrons from Z decay  
    vfs.addVetoOnThisFinalState(electronsFromZ);
    //addProjection(vfs, "JetFS");
    //jets
    D0ILConeJets jets(vfs);
    addProjection(jets, "Jets");
    //vertex
    PVertex vertex;
    addProjection(vertex, "PrimaryVertex");
    //jet isolation
    D0JetFromParticleIso jetiso(jets, electronsFromZ, new D0JetIsoEstimator(0.4), 20.);
    addProjection(jetiso, "JetIsolation");
  } 



  // Book histograms
  void D0_2008_S6879055::init() {
    // Use histogram auto-booking
    _crossSectionRatio = bookHistogram1D("crossSectionRatio", "$\\sigma(Z/\\gamma + >= n \\text{ jets}) / \\sigma(Z/\\gamma \\text{ inclusive})$", 5, -0.5, 4.5);
    //_crossSectionRatio =  bookHistogram1D(1,1,2, "#sigma(Z/#gamma + >=njets)/#sigma(Z/#gamma inclusive)");
    _crossSectionRatioNormToDataBin1 = bookHistogram1D("crossSectionRatioNormToDataBin1", "\\sigma(Z/\\gamma + >= n \\text{ jets}) / \\sigma(Z/\\gamma \\text{ inclusive})$", 5, -0.5, 4.5);
  }



  // Do the analysis 
  void D0_2008_S6879055::analyze(const Event & event) {
    getLog() << Log::DEBUG << "Starting analyzing" << endl;
    double weight = event.weight();
    getLog() << Log::DEBUG << "Event weight is " << weight << endl;

    //skip if the event is empty
    const FinalState & fs = applyProjection<FinalState>(event, "FS");
    if (fs.isEmpty()) {
      return;
    }

    //check that the primary vertex is within 60 cm in z from (0,0,0)
    const PVertex& vertex = applyProjection<PVertex>(event, "PrimaryVertex");
    getLog() << Log::DEBUG << "Primary vertex is at " << vertex.getPVPosition()/cm << " cm" << endl;
    // Skip event if the vertex position is not within 60 cm.
    if (fabs(vertex.getPVPosition().z())/cm > 60) {
      getLog() << Log::DEBUG << "Skipping event " << event.genEvent().event_number()
               << " because vertex position in z is " << vertex.getPVPosition().z()/cm << endl;
      return;
    }
    // Find the Z candidates
    const InvMassFinalState & invmassfs = applyProjection<InvMassFinalState>(event, "ElectronsFromZ");
    // If there is no Z candidate in the FinalState, skip the event
    if (invmassfs.isEmpty()) {
      getLog() << Log::DEBUG << "Skipping event " << event.genEvent().event_number()
               << " because no Z candidate is found " << endl;
      return;
    }
    // Now build the jets on a fs without the electrons from Z
    const D0ILConeJets & jetpro = applyProjection<D0ILConeJets>(event, "Jets");

    // Take the jets with pT > 20
    /// @todo Purge the evil!
    const list<FourMomentum>& jets = jetpro.getLorentzJets();
    list<FourMomentum> jets_aboveptmin;
    foreach (const FourMomentum& ijet, jets) {
      if (ijet.pT()/GeV > 20) {
        jets_aboveptmin.push_back(ijet);
      }
    }
    getLog() << Log::DEBUG << "Num jets above 20 GeV = " << jets_aboveptmin.size() << endl;

    // Check they are isolated from leptons
    const D0JetFromParticleIso & isoJet = applyProjection<D0JetFromParticleIso>(event, "JetIsolation");
    if (isoJet.getIsolatedParticles(1).size() != jets_aboveptmin.size()) {
      getLog() << Log::DEBUG << "Skipping event " << event.genEvent().event_number()
               << " because jet isolated from lepton size " << isoJet.getIsolatedParticles(0).size()
               << " is different from jet size " << jets_aboveptmin.size() << endl;
      return;
    }
    // Now take only the jets with |eta|<2.5
    list<FourMomentum> finaljet_list;
    foreach (const FourMomentum& ijet, jets_aboveptmin) {
      if (fabs(ijet.pseudorapidity()) < 2.5) {
        finaljet_list.push_back(ijet);
      }
    }
    getLog() << Log::DEBUG << "Num jets above 20 GeV and with |eta| < 2.5 = " << finaljet_list.size() << endl;

    // Now count....
    _events += weight;
    //_crossSectionRatio->fill(0, weight);
    if (finaljet_list.size() >= 1) {
      _crossSectionRatio->fill(1, weight);
      _crossSectionRatioNormToDataBin1->fill(1, weight);
    }
    if (finaljet_list.size() >= 2) {
      _crossSectionRatio->fill(2, weight);
      _crossSectionRatioNormToDataBin1->fill(2, weight);
    }
    if (finaljet_list.size() >= 3) {
      _crossSectionRatio->fill(3, weight);
      _crossSectionRatioNormToDataBin1->fill(3, weight);
    }
    if (finaljet_list.size() >= 4) {
      _crossSectionRatio->fill(4, weight);
      _crossSectionRatioNormToDataBin1->fill(4, weight);
    }
  }



  // Finalize
  void D0_2008_S6879055::finalize() {
    // Now divide by the inclusive result
    getLog() << Log::DEBUG << "Entries for bin 0 " << _crossSectionRatio->binEntries(0) << endl;
    if (_events != 0) {
      _crossSectionRatio->scale((double) 1 / _events);
      _crossSectionRatioNormToDataBin1->scale((double) 1 / _events);
      _crossSectionRatioNormToDataBin1->scale(0.1201 / _crossSectionRatio->binHeight(1));
    }
    getLog() << Log::DEBUG << "Finished!" << endl;
  }

}
