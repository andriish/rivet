// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {


  /// @brief Charged particle multiplicities and distributions
  class TASSO_1980_I143691 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TASSO_1980_I143691);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      ChargedFinalState cfs;
      declare(cfs, "CFS");
      declare(Thrust(cfs), "Thrust");
      

      // Book histograms
      _mult = bookCounter("/TMP/mult");
      unsigned int iloc(0);
      if(fuzzyEquals(sqrtS()/GeV, 13 , 1E-3))
	iloc = 1;
      else if(inRange(sqrtS(),16.99,23.01))
	iloc = 2;
      else if(inRange(sqrtS(),27.3,31.7))
	iloc = 3;
      else
	MSG_ERROR("Beam energy not supported!");

      _h_rap = bookHisto1D(iloc+1,1,1);
      _h_x   = bookHisto1D(iloc+4,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      // thrust
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      Vector3 axis=thrust.thrustAxis();
      foreach (const Particle& p, cfs.particles()) {
	const Vector3 mom3 = p.p3();
	double pp = mom3.mod();
	double xP = 2.*pp/sqrtS();
	_h_x->fill(xP,weight);
        const double mom = dot(axis, mom3);
	const double rap = 0.5 * log((p.E() + mom) /
				     (p.E() - mom));
	_h_rap->fill(fabs(rap),event.weight());
	_mult->fill(weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_rap, 1./sumOfWeights());
      scale(_h_x  , crossSection()*sqr(sqrtS())/sumOfWeights()/microbarn);  

      scale(_mult,1./sumOfWeights());
      
      Scatter2D temphisto(refData(1, 1, 1));
      Scatter2DPtr     mult = bookScatter2D(1, 1, 1);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
	const double x  = temphisto.point(b).x();
	pair<double,double> ex = temphisto.point(b).xErrs();
	pair<double,double> ex2 = ex;
	if(ex2.first ==0.) ex2. first=0.2;
	if(ex2.second==0.) ex2.second=0.2;
	if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
	  mult   ->addPoint(x, _mult->val(), ex, make_pair(_mult->err(), _mult->err()));
	}
	else {
	  mult   ->addPoint(x, 0., ex, make_pair(0.,.0));
	}
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_rap, _h_x;
    CounterPtr _mult;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TASSO_1980_I143691);


}
