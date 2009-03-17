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


  D0_2008_S6879055::D0_2008_S6879055() {
    setBeams(PROTON, ANTIPROTON);

    // Basic final state
    FinalState fs(-5.0, 5.0);
    addProjection(fs, "FS");

    // Leading electrons in tracking acceptance
    LeadingParticlesFinalState lpfs(fs, -1.1, 1.1, 25*GeV);
    lpfs.addParticleId(ELECTRON).addParticleId(POSITRON);
    addProjection(lpfs, "LeadingElectronsFS");

    // Invariant mass selection around Z pole
    InvMassFinalState electronsFromZ(lpfs, make_pair(ELECTRON, POSITRON), 75*GeV, 105*GeV);
    addProjection(electronsFromZ,"ElectronsFromZ");

    // Vetoed FS for jets
    VetoedFinalState vfs(fs);
    // Add particle/antiparticle vetoing
    vfs.vetoNeutrinos();
    // Veto the electrons from Z decay  
    vfs.addVetoOnThisFinalState(electronsFromZ);
    addProjection(vfs, "JetFS");

    // Jet finder
    D0ILConeJets jets(vfs, 0.5, 20.0*GeV);
    addProjection(jets, "Jets");

    // Vertex
    PVertex vertex;
    addProjection(vertex, "PrimaryVertex");
  } 



  // Book histograms
  void D0_2008_S6879055::init() {
    _crossSectionRatio = bookHistogram1D
      (1, 1, 1, "$\\sigma(Z/\\gamma + >= n \\text{ jets}) / \\sigma(Z/\\gamma \\text{ inclusive})$");
    _pTjet1 = bookHistogram1D
      (2, 1, 1, "$p_\\perp$ of 1st jet for $N_{\\mathrm{jet}} \\geq 1$ sample");
    _pTjet2 = bookHistogram1D
      (3, 1, 1, "$p_\\perp$ of 2nd jet for $N_{\\mathrm{jet}} \\geq 2$ sample");
    _pTjet3 = bookHistogram1D
      (4, 1, 1, "$p_\\perp$ of 3rd jet for $N_{\\mathrm{jet}} \\geq 3$ sample");
  }



  // Do the analysis 
  void D0_2008_S6879055::analyze(const Event& event) {
    const double weight = event.weight();

    // Skip if the event is empty
    const FinalState& fs = applyProjection<FinalState>(event, "FS");
    if (fs.isEmpty()) {
      vetoEvent(event);
    }

    // Check that the primary vertex is within 60 cm in z from (0,0,0)
    const PVertex& vertex = applyProjection<PVertex>(event, "PrimaryVertex");
    getLog() << Log::DEBUG << "Primary vertex is at " << vertex.position()/cm << " cm" << endl;
    if (fabs(vertex.position().z())/cm > 60) {
      getLog() << Log::DEBUG << "Vertex z-position " << vertex.position().z()/cm << " is outside cuts" << endl;
      vetoEvent(event);
    }

    // Find the Z candidates
    const InvMassFinalState& invmassfs = applyProjection<InvMassFinalState>(event, "ElectronsFromZ");
    // If there is no Z candidate in the FinalState, skip the event
    if (invmassfs.particles().size()!=2) {
      getLog() << Log::DEBUG << "No Z candidate found" << endl;
      vetoEvent(event);
    }

    // Now build the jets on a FS without the electrons from Z
    /// @todo Purge the evil!
    const D0ILConeJets& jetpro = applyProjection<D0ILConeJets>(event, "Jets");

    // additional cuts on jets: |eta| < 2.5 and dR(j,leading electron) > 0.4
    const list<FourMomentum>& jets = jetpro.getLorentzJets();
    vector<FourMomentum> finaljet_list;
    foreach (const FourMomentum& ijet, jets) {
      FourMomentum e0=invmassfs.particles()[0].momentum();
      FourMomentum e1=invmassfs.particles()[0].momentum();
      if (fabs(ijet.pseudorapidity()) < 2.5 &&
          deltaR(e0.pseudorapidity(), e0.azimuthalAngle(),
                 ijet.pseudorapidity(), ijet.azimuthalAngle()) > 0.4 &&
          deltaR(e1.pseudorapidity(), e1.azimuthalAngle(),
                 ijet.pseudorapidity(), ijet.azimuthalAngle()) > 0.4) {
        finaljet_list.push_back(ijet);
      }
    }
    getLog() << Log::DEBUG << "Num jets passing = " << finaljet_list.size() << endl;

    // For normalisation of crossSection data
    _crossSectionRatio->fill(0, weight);

    // Fill jet pT and multiplicities
    if (finaljet_list.size() >= 1) {
      _crossSectionRatio->fill(1, weight);
      _pTjet1->fill(finaljet_list[0].pT(), weight);
    }
    if (finaljet_list.size() >= 2) {
      _crossSectionRatio->fill(2, weight);
      _pTjet2->fill(finaljet_list[1].pT(), weight);
    }
    if (finaljet_list.size() >= 3) {
      _crossSectionRatio->fill(3, weight);
      _pTjet3->fill(finaljet_list[2].pT(), weight);
    }
    if (finaljet_list.size() >= 4) {
      _crossSectionRatio->fill(4, weight);
    }
  }



  // Finalize
  void D0_2008_S6879055::finalize() {
    // Now divide by the inclusive result
    _crossSectionRatio->scale(1.0/_crossSectionRatio->binHeight(0));

    // Normalise jet pT's to integral of data
    normalize(_pTjet1, 10439.0);
    normalize(_pTjet2, 1461.5);
    normalize(_pTjet3, 217.0);
  }


}
