// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Forward photon production cross-section at 13 TeV
  class LHCF_2018_I1518782 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(LHCF_2018_I1518782);

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

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      Particles fs_particles = apply<FinalState> (event, "FS").particles();
      for (const Particle& p : fs_particles) {

        // Select photons above threshold
        if (p.abspid() != PID::PHOTON) continue;
        if (p.E() < 200*GeV) continue;

        // Double analysis efficiency with a two-sided LHCf
        const double eta = abs(p.eta());
        const double energy = p.E()/GeV;

        // Fill histograms
        if ( eta > 10.94 ) {
          _h_n_en_eta1->fill(energy);
        }
        else if (eta > 8.81 && eta < 8.99) {
          _h_n_en_eta2->fill(energy);
        }

      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // Norm to cross-section
      scale(_h_n_en_eta1, crossSection()/millibarn/sumOfWeights()/2.);
      scale(_h_n_en_eta2, crossSection()/millibarn/sumOfWeights()/2.);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_n_en_eta1, _h_n_en_eta2;
    //@}

  };


  DECLARE_RIVET_PLUGIN(LHCF_2018_I1518782);

}
