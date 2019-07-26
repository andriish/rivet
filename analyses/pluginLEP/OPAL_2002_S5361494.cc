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
  class OPAL_2002_S5361494 : public Analysis {
  public:

    /// Constructor
    OPAL_2002_S5361494()
      : Analysis("OPAL_2002_S5361494"),
        _weightLight(0.),_weightCharm(0.),_weightBottom(0.)
    {}

    /// @name Analysis methods
    //@{


    void init() {
      // Projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "CFS");
      declare(InitialQuarks(), "IQF");
      _cLight  = bookCounter("/TMP/CLIGHT" );
      _cCharm  = bookCounter("/TMP/CCHARM" );
      _cBottom = bookCounter("/TMP/CBOTTOM");
    }


    void analyze(const Event& event) {
      const double weight = event.weight();
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
      const size_t numParticles = cfs.particles().size();
      switch (flavour) {
      case 1: case 2: case 3:
	_weightLight  += weight;
        _cLight->fill(  numParticles * weight );
        break;
      case 4:
        _weightCharm  += weight;
	_cCharm->fill(  numParticles * weight );
        break;
      case 5:
        _weightBottom += weight;
        _cBottom->fill( numParticles * weight );
        break;
      }

    }


    void finalize() {
      // calculate the averages and diffs
      if(_weightLight!=0)  scale( _cLight, 1./_weightLight);
      if(_weightCharm !=0) scale( _cCharm, 1./_weightCharm);
      if(_weightBottom!=0) scale(_cBottom,1./_weightBottom);
      Counter _cDiff = *_cBottom - *_cLight;
      // fill the histograms
      for(unsigned int ix=1;ix<5;++ix) {
	double val(0.), err(0.0);
	if(ix==1) {
	  val = _cBottom->val();
	  err = _cBottom->err();
	}
	else if(ix==2) {
	  val = _cCharm->val();
	  err = _cCharm->err();
	}
	else if(ix==3) {
	  val = _cLight->val();
	  err = _cLight->err();
	}
	else if(ix==4) {
	  val = _cDiff.val();
	  err = _cDiff.err();
	}
	Scatter2D temphisto(refData(1, 1, ix));
	Scatter2DPtr  mult = bookScatter2D(1, 1, ix);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
	    mult->addPoint(x, val, ex, make_pair(err,err));
	  }
	  else {
	    mult->addPoint(x, 0., ex, make_pair(0.,.0));
	  }
	}
      }
    }

    //@}


  private:

    /// @name Multiplicities
    //@{
    CounterPtr _cLight;
    CounterPtr _cCharm;
    CounterPtr _cBottom;
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
