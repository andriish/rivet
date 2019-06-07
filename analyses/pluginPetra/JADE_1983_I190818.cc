// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Average multiplcity at a range of energies
  class JADE_1983_I190818 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(JADE_1983_I190818);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      if( !(fuzzyEquals(sqrtS()/GeV,12.0) ||
	    fuzzyEquals(sqrtS()/GeV,30.0) ||
	    fuzzyEquals(sqrtS()/GeV,35.0) )) {
        MSG_WARNING("CoM energy of events sqrt(s) = " << sqrtS()/GeV
                    << " doesn't match any available analysis energy .");
      }
      book(_counter, "/TMP/MULT");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      MSG_DEBUG("Total charged multiplicity = " << cfs.size());
      _counter->fill(cfs.size());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_counter,1./sumOfWeights());

      double val = _counter->val();
      double err = _counter->err();
      
      Scatter2D temphisto(refData(1, 1, 1));
      Scatter2DPtr mult;
      book(mult, 1, 1, 1);
      
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
	  mult->addPoint(x,  0., ex,   make_pair(0.,.0));
	}
      }
    }
    //@}

  private:

    // Histogram
    CounterPtr _counter;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(JADE_1983_I190818);

}
