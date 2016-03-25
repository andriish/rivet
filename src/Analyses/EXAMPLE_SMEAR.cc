// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/SmearedJets.hh"

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

      // SmearedJets sj2(fj, JET_SMEAR_IDENTITY,
      //                 [](const Jet& j){ return j.bTagged() ? 0.7*(1 - exp(-j.pT()/(10*GeV))) : 0; } );
      // addProjection(sj2, "Jets2");

      // SmearedJets sj3(fj,
      //                 [](const Jet& j) { return j; },
      //                 [](const Jet& j){ return j.bTagged() ? 0.7*(1 - exp(-j.pT()/(10*GeV))) : 0; } );
      // addProjection(sj3, "Jets3");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const Jets jets0 = applyProjection<JetAlg>(event, "Jets0").jets(Cuts::pT > 10*GeV);
      const Jets jets1 = applyProjection<JetAlg>(event, "Jets1").jets(Cuts::pT > 10*GeV);
      const Jets jets2 = applyProjection<JetAlg>(event, "Jets2").jets(Cuts::pT > 10*GeV);
      const Jets jets3 = applyProjection<JetAlg>(event, "Jets3").jets(Cuts::pT > 10*GeV);
      MSG_INFO("Numbers of jets = " << jets0.size() << " true; "
               << jets1.size() << ", " << jets2.size() << ", " << jets3.size());
    }


    // /// Normalise histograms etc., after the run
    // void finalize() {
    //   /// @todo Normalise, scale and otherwise manipulate histograms here
    //   // scale(_h_YYYY, crossSection()/sumOfWeights()); // norm to cross section
    //   // normalize(_h_YYYY); // normalize to unity
    // }

    //@}


  private:

    // /// @name Histograms
    // //@{
    // Profile1DPtr _h_XXXX;
    // Histo1DPtr _h_YYYY;
    // //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(EXAMPLE_SMEAR);


}
