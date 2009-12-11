// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Measurement of the ratio sigma(Z/gamma* + n jets)/sigma(Z/gamma*)
  class D0_2008_S6879055 : public Analysis {
  public:
 
    /// Default constructor.
    D0_2008_S6879055() : Analysis("D0_2008_S6879055")
    {
      setBeams(PROTON, ANTIPROTON);
    }


    /// @name Analysis methods
    //@{
 
    // Book histograms
    void init() {
      // Basic final state
      FinalState fs(-5.0, 5.0);
      addProjection(fs, "FS");
   
      // Leading electrons in tracking acceptance
      LeadingParticlesFinalState lpfs(FinalState(-1.1, 1.1, 25*GeV));
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
      FastJets jets(vfs, FastJets::D0ILCONE, 0.5, 20.0*GeV);
      addProjection(jets, "Jets");
   
      // Vertex
      PVertex vertex;
      addProjection(vertex, "PrimaryVertex");

      _crossSectionRatio = bookHistogram1D(1, 1, 1);
      _pTjet1 = bookHistogram1D(2, 1, 1);
      _pTjet2 = bookHistogram1D(3, 1, 1);
      _pTjet3 = bookHistogram1D(4, 1, 1);
    }
 
 
 
    /// Do the analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
   
      // Skip if the event is empty
      const FinalState& fs = applyProjection<FinalState>(event, "FS");
      if (fs.empty()) {
        vetoEvent;
      }
   
      // Check that the primary vertex is within 60 cm in z from (0,0,0)
      const PVertex& vertex = applyProjection<PVertex>(event, "PrimaryVertex");
      getLog() << Log::DEBUG << "Primary vertex is at " << vertex.position()/cm << " cm" << endl;
      if (fabs(vertex.position().z())/cm > 60) {
        getLog() << Log::DEBUG << "Vertex z-position " << vertex.position().z()/cm << " is outside cuts" << endl;
        vetoEvent;
      }
   
      // Find the Z candidates
      const InvMassFinalState& invmassfs = applyProjection<InvMassFinalState>(event, "ElectronsFromZ");
      // If there is no Z candidate in the FinalState, skip the event
      if (invmassfs.particles().size() != 2) {
        getLog() << Log::DEBUG << "No Z candidate found" << endl;
        vetoEvent;
      }
   
      // Now build the list of jets on a FS without the electrons from Z
      // Additional cuts on jets: |eta| < 2.5 and dR(j,leading electron) > 0.4
      const JetAlg& jetpro = applyProjection<JetAlg>(event, "Jets");
      const Jets jets = jetpro.jetsByPt(20.0*GeV);
      vector<FourMomentum> finaljet_list;
      foreach (const Jet& j, jets) {
        const double jeta = j.momentum().eta();
        const double jphi = j.momentum().phi();
        if (fabs(jeta) > 2.5) continue;
     
        FourMomentum e0 = invmassfs.particles()[0].momentum();
        FourMomentum e1 = invmassfs.particles()[1].momentum();
        const double e0eta = e0.pseudorapidity();
        const double e0phi = e0.azimuthalAngle();
        if (deltaR(e0eta, e0phi, jeta, jphi) < 0.4) continue;
     
        const double e1eta = e1.pseudorapidity();
        const double e1phi = e1.azimuthalAngle();
        if (deltaR(e1eta, e1phi, jeta, jphi) < 0.4) continue;
     
        // If we pass all cuts...
        finaljet_list.push_back(j.momentum());
      }
      getLog() << Log::DEBUG << "Num jets passing = " << finaljet_list.size() << endl;
   
      // For normalisation of crossSection data (includes events with no jets passing cuts)
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
 
 
 
    /// Finalize
    void finalize() {
      // Now divide by the inclusive result
      _crossSectionRatio->scale(1.0/_crossSectionRatio->binHeight(0));
   
      // Normalise jet pT's to integral of data
      // NB. There is no other way to do this, because these quantities are not
      // detector-corrected
      normalize(_pTjet1, 10439.0);
      normalize(_pTjet2, 1461.5);
      normalize(_pTjet3, 217.0);
    }
 
    //@}


  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D * _crossSectionRatio;
    AIDA::IHistogram1D * _pTjet1;
    AIDA::IHistogram1D * _pTjet2;
    AIDA::IHistogram1D * _pTjet3;
    //@}

  };

 
 
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<D0_2008_S6879055> plugin_D0_2008_S6879055;

}
