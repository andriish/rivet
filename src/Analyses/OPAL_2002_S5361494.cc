// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
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
  class OPAL_2002_S5361494 : public Analysis {
  public:

    /// Constructor
    OPAL_2002_S5361494()
      : Analysis("OPAL_2002_S5361494"),
	_weightedTotalChargedPartNumLight(0.),
	_weightedTotalChargedPartNumCharm(0.),
	_weightedTotalChargedPartNumBottom(0.),
	_weightLight(0.),_weightCharm(0.),_weightBottom(0.)
    {}

    /// @name Analysis methods
    //@{


    void init() {
      // Projections
      addProjection(Beam(), "Beams");
      addProjection(ChargedFinalState(), "CFS");
      addProjection(InitialQuarks(), "IQF");

    }


    void analyze(const Event& event) {
      const double weight = event.weight();
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      const FinalState& cfs = applyProjection<FinalState>(event, "CFS");
      if (cfs.size() < 2) vetoEvent;


      int flavour = 0;
      const InitialQuarks& iqf = applyProjection<InitialQuarks>(event, "IQF");
      
      // If we only have two quarks (qqbar), just take the flavour.
      // If we have more than two quarks, look for the highest energetic q-qbar pair.
      if (iqf.particles().size() == 2) {
        flavour = abs( iqf.particles().front().pdgId() );
      }
      else {
        map<int, double> quarkmap;
        foreach (const Particle& p, iqf.particles()) {
          if (quarkmap[p.pdgId()] < p.momentum().E()) {
            quarkmap[p.pdgId()] = p.momentum().E();
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
	_weightLight  += weight;
        _weightedTotalChargedPartNumLight  += numParticles * weight;
        break;
      case 4:
	_weightCharm  += weight;
        _weightedTotalChargedPartNumCharm  += numParticles * weight;
        break;
      case 5:
	_weightBottom += weight;
        _weightedTotalChargedPartNumBottom += numParticles * weight;
        break;
      }

    }


    void finalize() {
      // bottom
      const double avgNumPartsBottom = _weightedTotalChargedPartNumBottom / _weightBottom;
      AIDA::IDataPointSet * multB = bookDataPointSet(1, 1, 1);
      for (int i = 0; i < multB->size(); ++i) {
        if (fuzzyEquals(sqrtS(), multB->point(i)->coordinate(0)->value(), 0.01)) {
          multB->point(i)->coordinate(1)->setValue(avgNumPartsBottom);
        }
      }
      // charm
      const double avgNumPartsCharm = _weightedTotalChargedPartNumCharm / _weightCharm;
      AIDA::IDataPointSet * multC = bookDataPointSet(1, 1, 2);
      for (int i = 0; i < multC->size(); ++i) {
        if (fuzzyEquals(sqrtS(), multC->point(i)->coordinate(0)->value(), 0.01)) {
          multC->point(i)->coordinate(1)->setValue(avgNumPartsCharm);
        }
      }
      // light
      const double avgNumPartsLight = _weightedTotalChargedPartNumLight / _weightLight;
      AIDA::IDataPointSet * multL = bookDataPointSet(1, 1, 3);
      for (int i = 0; i < multL->size(); ++i) {
        if (fuzzyEquals(sqrtS(), multL->point(i)->coordinate(0)->value(), 0.01)) {
          multL->point(i)->coordinate(1)->setValue(avgNumPartsLight);
        }
      }
      // bottom-light
      AIDA::IDataPointSet * multD = bookDataPointSet(1, 1, 4);
      for (int i = 0; i < multD->size(); ++i) {
        if (fuzzyEquals(sqrtS(), multD->point(i)->coordinate(0)->value(), 0.01)) {
          multD->point(i)->coordinate(1)->setValue(avgNumPartsBottom-avgNumPartsLight);
        }
      }
    }
    //@}


  private:

    /// @name Multiplicities
    //@{
    double _weightedTotalChargedPartNumLight;
    double _weightedTotalChargedPartNumCharm;
    double _weightedTotalChargedPartNumBottom;
    //@}

    /// @name Weights
    //@{
    double _weightLight;
    double _weightCharm;
    double _weightBottom;
    //@}
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_2002_S5361494);

}
