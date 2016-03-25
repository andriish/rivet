// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/SmearedJets.hh"

namespace Rivet {


  class SMEAR : public Analysis {
  public:

    /// Constructor
    SMEAR()
      : Analysis("SMEAR")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      FastJets fj(FinalState(Cuts::abseta < 5), FastJets::ANTIKT, 0.4);
      SmearedJets sj(fj, JET_EFF_ONE, JET_SMEAR_IDENTITY);
      // SmearedJets sj(fj, [](const Jet& j) -> double { return 1 - exp(-j.pT()/(10*GeV)); });
      // SmearedJets sj(fj,
      //                [](const Jet& j) { return 1 - exp(-j.pT()/(10*GeV)); },
      //                [](const Jet& j) { return j; });
      addProjection(sj, "Jets");
      SmearedJets sj2(fj, JET_EFF_ZERO, JET_SMEAR_IDENTITY);
      addProjection(sj, "Jets0");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const Jets jets = applyProjection<JetAlg>(event, "Jets").jets(Cuts::pT > 10*GeV);
      MSG_INFO("Number of jets = " << jets.size());
      const Jets jets0 = applyProjection<JetAlg>(event, "Jets0").jets(Cuts::pT > 10*GeV);
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
  DECLARE_RIVET_PLUGIN(SMEAR);


}
