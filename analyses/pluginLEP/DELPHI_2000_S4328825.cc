// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/InitialQuarks.hh"
#include <cmath>

namespace Rivet {


  /// @brief OPAL multiplicities at various energies
  /// @author Peter Richardson
  class DELPHI_2000_S4328825 : public Analysis {
  public:

    /// Constructor
    DELPHI_2000_S4328825()
      : Analysis("DELPHI_2000_S4328825")
    {}

    /// @name Analysis methods
    //@{


    void init() {
      // Projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "CFS");
      declare(InitialQuarks(), "IQF");
      book(_weightedTotalChargedPartNumLight,"weight_totalch_light");
      book(_weightedTotalChargedPartNumCharm,"weight_totalch_charm");
      book(_weightedTotalChargedPartNumBottom,"weight_totalch_bottom");
      book(_weightLight,"weight_light");
      book(_weightCharm,"weight_charm");
      book(_weightBottom,"weight_bottom");

      book(h_bottom, 1, 1, 1);
      book(h_charm, 1, 1, 2);
      book(h_light, 1, 1, 3);
      book(h_diff, 1, 1, 4);  // bottom minus light

    }


    void analyze(const Event& event) {
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      if (cfs.size() < 2) vetoEvent;


      int flavour = 0;
      const InitialQuarks& iqf = apply<InitialQuarks>(event, "IQF");

      // If we only have two quarks (qqbar), just take the flavour.
      // If we have more than two quarks, look for the highest energetic q-qbar pair.
      if (iqf.particles().size() == 2) {
        flavour = iqf.particles().front().abspid();
      }
      else {
        map<int, double> quarkmap;
        for (const Particle& p : iqf.particles()) {
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
      const size_t numParticles = cfs.particles().size();
      switch (flavour) {
      case 1: case 2: case 3:
        _weightLight->fill();
        _weightedTotalChargedPartNumLight ->fill(numParticles);
        break;
      case 4:
        _weightCharm->fill();
        _weightedTotalChargedPartNumCharm ->fill(numParticles);
        break;
      case 5:
        _weightBottom->fill();
        _weightedTotalChargedPartNumBottom->fill(numParticles);
        break;
      }

    }


    void finalize() {
      Histo1D temphisto(refData(1, 1, 1));

      const double avgNumPartsBottom = dbl(*_weightedTotalChargedPartNumBottom / *_weightBottom);
      const double avgNumPartsCharm  = dbl(*_weightedTotalChargedPartNumCharm  / *_weightCharm);
      const double avgNumPartsLight  = dbl(*_weightedTotalChargedPartNumLight  / *_weightLight);
      for (size_t b = 0; b < temphisto.numBins(); b++) {
        const double x  = temphisto.bin(b).xMid();
        const double ex = temphisto.bin(b).xWidth()/2.;
        if (inRange(sqrtS()/GeV, x-ex, x+ex)) {
          // @TODO: Fix y-error:
          h_bottom->addPoint(x, avgNumPartsBottom, ex, 0.);
          h_charm->addPoint(x, avgNumPartsCharm, ex, 0.);
          h_light->addPoint(x, avgNumPartsLight, ex, 0.);
          h_diff->addPoint(x, avgNumPartsBottom-avgNumPartsLight, ex, 0.);
        }
      }
    }

    //@}


  private:
    Scatter2DPtr h_bottom, h_charm, h_light, h_diff;


    /// @name Multiplicities
    //@{
    CounterPtr _weightedTotalChargedPartNumLight;
    CounterPtr _weightedTotalChargedPartNumCharm;
    CounterPtr _weightedTotalChargedPartNumBottom;
    //@}

    /// @name Weights
    //@{
    CounterPtr _weightLight;
    CounterPtr _weightCharm;
    CounterPtr _weightBottom;
    //@}
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(DELPHI_2000_S4328825);

}
