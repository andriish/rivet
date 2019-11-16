// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  /// @brief charged particle spectra
  class DELPHI_1997_I428178 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(DELPHI_1997_I428178);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(InitialQuarks(), "IQF");


      // Book histograms
      _h_bottom = bookHisto1D(1, 1, 1);
      _h_charm  = bookHisto1D(1, 1, 2);
      _h_light  = bookHisto1D(1, 1, 3);

      _wBottom = _wCharm = _wLight = 0.;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // Get event weight for histo filling
      const double weight = event.weight();

      int flavour = 0;
      const InitialQuarks& iqf = apply<InitialQuarks>(event, "IQF");

      // If we only have two quarks (qqbar), just take the flavour.
      // If we have more than two quarks, look for the highest energetic q-qbar pair.
      if (iqf.particles().size() == 2) {
        flavour = iqf.particles().front().abspid();
      }
      else {
        map<int, double> quarkmap;
        foreach (const Particle& p, iqf.particles()) {
          if (quarkmap[p.pid()] < p.E()) {
            quarkmap[p.pid()] = p.E();
          }
        }
        double maxenergy = 0.;
        for (int i = 1; i <= 5; ++i) {
          if (quarkmap[i]+quarkmap[-i] > maxenergy) {
            flavour = i;
          }
        }
      }
      if     (flavour==5) _wBottom += weight;
      else if(flavour==4) _wCharm  += weight;
      else                _wLight  += weight;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      foreach (const Particle& p, fs.particles()) {
        double xp = p.p3().mod()/meanBeamMom;
	if     (flavour==5) _h_bottom->fill(xp,weight);
	else if(flavour==4) _h_charm ->fill(xp,weight);
	else                _h_light ->fill(xp,weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_bottom, 1./_wBottom);
      scale(_h_charm , 1./_wCharm);
      scale(_h_light , 1./_wLight);

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_bottom, _h_charm, _h_light;
    double _wBottom, _wCharm, _wLight;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(DELPHI_1997_I428178);


}
