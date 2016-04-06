// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/SmearedJets.hh"
#include "Rivet/Projections/SmearedParticles.hh"

namespace Rivet {


  class EXAMPLE_SMEAR : public Analysis {
  public:

    /// Constructor
    EXAMPLE_SMEAR()
      : Analysis("EXAMPLE_SMEAR")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      FastJets fj(FinalState(Cuts::abseta < 5), FastJets::ANTIKT, 0.4);
      addProjection(fj, "Jets0");

      SmearedJets sj1(fj, JET_SMEAR_IDENTITY);
      addProjection(sj1, "Jets1");

      SmearedJets sj2(fj, JET_SMEAR_ATLAS_RUN1,
                      [](const Jet& j){ return j.bTagged() ? 0.7*(1 - exp(-j.pT()/(10*GeV))) : 0; } );
      addProjection(sj2, "Jets2");

      SmearedJets sj3(fj,
                      [](const Jet& j) { return j; },
                      [](const Jet& j){ return j.bTagged() ? 0.7*(1 - exp(-j.pT()/(10*GeV))) : 0; },
                      JET_CTAG_PERFECT,
                      JET_EFF_ZERO);
      addProjection(sj3, "Jets3");

      IdentifiedFinalState truthelectrons(Cuts::abseta < 2.5 && Cuts::pT > 10*GeV, {{PID::ELECTRON, PID::POSITRON}});
      addProjection(truthelectrons, "Electrons0");
      SmearedParticles recoelectrons(truthelectrons, ELECTRON_EFF_ATLAS_RUN1);
      addProjection(recoelectrons, "Electrons1");

      _h_njtrue = bookHisto1D("njets_true", 10, -0.5, 9.5);
      _h_njreco = bookHisto1D("njets_reco", 10, -0.5, 9.5);
      _h_netrue = bookHisto1D("nelec_true", 5, -0.5, 4.5);
      _h_nereco = bookHisto1D("nelec_reco", 5, -0.5, 4.5);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const Jets jets0 = applyProjection<JetAlg>(event, "Jets0").jets(Cuts::pT > 10*GeV);
      const Jets jets1 = applyProjection<JetAlg>(event, "Jets1").jets(Cuts::pT > 10*GeV);
      const Jets jets2 = applyProjection<JetAlg>(event, "Jets2").jets(Cuts::pT > 10*GeV);
      const Jets jets3 = applyProjection<JetAlg>(event, "Jets3").jets(Cuts::pT > 10*GeV);
      MSG_INFO("Numbers of jets = " << jets0.size() << " true; "
               << jets1.size() << ", " << jets2.size() << ", " << jets3.size());
      _h_njtrue->fill(jets0.size());
      _h_njreco->fill(jets1.size());

      const Particles& elecs0 = applyProjection<ParticleFinder>(event, "Electrons0").particles();
      const Particles& elecs1 = applyProjection<ParticleFinder>(event, "Electrons1").particles();
      MSG_INFO("Numbers of electrons = " << elecs0.size() << " true; " << elecs1.size() << " reco");
      _h_netrue->fill(elecs0.size());
      _h_nereco->fill(elecs1.size());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_njtrue);
      normalize(_h_njreco);
      normalize(_h_netrue);
      normalize(_h_nereco);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_njtrue, _h_njreco, _h_netrue, _h_nereco;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(EXAMPLE_SMEAR);


}
