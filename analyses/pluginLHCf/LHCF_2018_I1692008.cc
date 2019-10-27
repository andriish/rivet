// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {

  /// @brief forward neutron production cross-section at 13 TeV
  class LHCF_2018_I1692008 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(LHCF_2018_I1692008);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beam");
      declare(FinalState(), "FS");

      // Book histograms
      book(_h_n_en_eta1, 1, 1, 1);
      book(_h_n_en_eta2, 2, 1, 1);
      book(_h_n_en_eta3, 3, 1, 1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const FinalState &fs = apply<FinalState> (event, "FS");
      Particles fs_particles = fs.particles();

      for (Particle& p : fs_particles) {

        // select photons above threshold
        if (p.abspid() != 2112) continue; //Select neutrons and antineutrons only
        if (p.E()/GeV < 500.)   continue; //Select only energy above 500 GeV
        //if(p.pT()/GeV == 0.0) continue; //Artificially remove QGSJet II-04 wrong events

        // Double analysis efficiency with a two-sided LHCf
        const double eta = abs(p.eta());
        const double energy = p.E()/GeV;

        if ( eta > 10.758267 ) {
          _h_n_en_eta1->fill( energy );
        } 
        else if (eta > 8.994669 && eta < 9.217812) {
          _h_n_en_eta2->fill( energy );
        } 
        else if (eta > 8.812347 && eta < 8.994669) {
          _h_n_en_eta3->fill( energy );
        }
      }

    }

    /// Normalise histograms etc., after the run
    void finalize() {
      //Scale considering the LHCf Arm2 side
      scale(_h_n_en_eta1, crossSection()/millibarn/sumOfWeights()/2.); // norm to cross section
      scale(_h_n_en_eta2, crossSection()/millibarn/sumOfWeights()/2.); // norm to cross section
      scale(_h_n_en_eta3, crossSection()/millibarn/sumOfWeights()/2.); // norm to cross section
    }
    //@}

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_n_en_eta1;
    Histo1DPtr _h_n_en_eta2;
    Histo1DPtr _h_n_en_eta3;
    //@}
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(LHCF_2018_I1692008);

}
