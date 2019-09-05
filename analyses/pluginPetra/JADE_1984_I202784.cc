// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief EEC at 14, 22 and 34 GeV
  class JADE_1984_I202784 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(JADE_1984_I202784);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      unsigned int iloc(0);
      if(fuzzyEquals(sqrtS()/GeV, 14., 1e-3)) {
	iloc=1;
      }
      else if(fuzzyEquals(sqrtS()/GeV, 22., 1e-3)) {
	iloc=2;
      }
      else if (fuzzyEquals(sqrtS()/GeV, 34., 1e-3)) {
	iloc=3;
      }
      else
	MSG_ERROR("Beam energy not supported!");
      _histEEC  = bookHisto1D(1, 1, iloc);
      _histAEEC = bookHisto1D(2, 1, iloc);
      _weightSum =0.;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if ( fs.particles().size() < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");
      const double weight = event.weight();
      _weightSum += weight;

      double Evis = 0.0;
      foreach (const Particle& p, fs.particles()) {
        Evis += p.E();
      }
      double Evis2 = sqr(Evis);
      // (A)EEC
      // Need iterators since second loop starts at current outer loop iterator, i.e. no "foreach" here!
      for (Particles::const_iterator p_i = fs.particles().begin(); p_i != fs.particles().end(); ++p_i) {
        for (Particles::const_iterator p_j = p_i; p_j != fs.particles().end(); ++p_j) {
          const Vector3 mom3_i = p_i->momentum().p3();
          const Vector3 mom3_j = p_j->momentum().p3();
          const double energy_i = p_i->momentum().E();
          const double energy_j = p_j->momentum().E();
          const double thetaij = mom3_i.unit().angle(mom3_j.unit())/M_PI*180.;
          double eec = (energy_i*energy_j) / Evis2;
	  if(p_i != p_j) eec *= 2.;
	  _histEEC ->fill(thetaij,  eec*weight);
          if (thetaij < 90.) {
            _histAEEC->fill(thetaij, -eec*weight);
	  }
          else {
            _histAEEC  ->fill(180.-thetaij, eec*weight);
	  }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_histEEC   , 180.0/M_PI/_weightSum*1000.);
      scale(_histAEEC  , 180.0/M_PI/_weightSum*1000.);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _histEEC, _histAEEC;
    double _weightSum;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(JADE_1984_I202784);


}
