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
    _crossSectionRatio = bookHistogram1D("crossSectionRatio", "#sigma(Z/#gamma + >=njets)/#sigma(Z/#gamma inclusive)", 5, -0.5, 4.5);
    //_crossSectionRatio =  bookHistogram1D(1,1,2, "#sigma(Z/#gamma + >=njets)/#sigma(Z/#gamma inclusive)");
    _crossSectionRatioNormToDataBin1 = bookHistogram1D("crossSectionRatioNormToDataBin1", "#sigma(Z/#gamma + >=njets)/#sigma(Z/#gamma inclusive)", 5, -0.5, 4.5);
  }

  // Do the analysis 
  void D0_2008_S6879055::analyze(const Event & event) {
    Log & log = getLog();
    log << Log::DEBUG << "Starting analyzing" << endl;
    double weight = event.weight();
    log << Log::DEBUG << "Event weight is " << weight << endl;

    //skip if the event is empty
    const FinalState & fs = applyProjection<FinalState>(event, "FS");
    if (fs.isEmpty())
      return;

    //check that the primary vertex is within 60 cm in z from (0,0,0)
    const PVertex & vertex = applyProjection<PVertex>(event, "PrimaryVertex");
    log << Log::DEBUG << "Primary vertex is at " << vertex.getPVPosition() << endl;
    // skip event if the vertex position is not within 60 cm.
    if (fabs(vertex.getPVPosition().z()) > 60.) {
      log << Log::DEBUG << "Skipping event " << event.genEvent().event_number()
        << " because vertex position in z is " << vertex.getPVPosition().z() << endl;
      return;
    }
    //find the Z candidates
    const InvMassFinalState & invmassfs = applyProjection<InvMassFinalState>(event, "ElectronsFromZ");
    //if there is no Z candidate in the FinalState skip event,
    if (invmassfs.isEmpty()) {
      log << Log::DEBUG << "Skipping event " << event.genEvent().event_number()
        << " because no Z candidate is found " << endl;
      return;
    }
    //now build the jets on a fs without the electrons from Z
    const D0ILConeJets & jetpro = applyProjection<D0ILConeJets>(event, "Jets");

    //take the jets with pt > 20
    const list < FourMomentum > &jets = jetpro.getLorentzJets();
    list < FourMomentum > jets_aboveptmin;
    list < FourMomentum >::const_iterator jetbegin = jets.begin();
    list < FourMomentum >::const_iterator jetend = jets.end();
    for (list < FourMomentum >::const_iterator ijet = jetbegin; ijet != jetend; ++ijet) {
      if (ijet->pT() > 20)
        jets_aboveptmin.push_back(*ijet);
    }
    log << Log::DEBUG << "jets above 20 GeV size = " << jets_aboveptmin.size() << endl;

    //check they are isolated form leptons
    const D0JetFromParticleIso & isoJet = applyProjection<D0JetFromParticleIso>(event, "JetIsolation");
    if (isoJet.getIsolatedParticles(1).size() != jets_aboveptmin.size()) {
      log << Log::DEBUG << "skipping event " << event.genEvent().event_number()
        << " because Jet Isolated from lepton size " << isoJet.getIsolatedParticles(0).size()
        << " is different from jet size " << jets_aboveptmin.size() << endl;
      return;
    }
    //now take only the jets with |eta|<2.5
    list < FourMomentum > finaljet_list;
    for (list < FourMomentum >::const_iterator ijet = jets_aboveptmin.begin(); ijet != jets_aboveptmin.end(); ++ijet) {
      if (fabs(ijet->pseudorapidity()) < 2.5)
        finaljet_list.push_back(*ijet);
      //if (fabs(ijet->pseudorapidity()) < 5) finaljet_list.push_back(*ijet);
    }
    log << Log::DEBUG << "jets above 20 GeV and with |eta| < 2.5  size = " << finaljet_list.size() << endl;

    //now count....
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

    // Finished
    log << Log::DEBUG << "Finished analyzing" << endl;
  }

  // Finalize
  void D0_2008_S6879055::finalize() {
    Log & log = getLog();

    //now divide by the inclusive result
    log << Log::DEBUG << "Entries for bin 0 " << _crossSectionRatio->binEntries(0) << endl;
    //if (_crossSectionRatio->binEntries(0) != 0){
    if (_events != 0) {
      _crossSectionRatio->scale((double) 1 / _events);
      _crossSectionRatioNormToDataBin1->scale((double) 1 / _events);
      _crossSectionRatioNormToDataBin1->scale(0.1201 / _crossSectionRatio->binHeight(1));
    }
    log << Log::DEBUG << "Finished!" << endl;

  }

}
