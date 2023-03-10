// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class DM2_1991_I318558 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DM2_1991_I318558);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      // Book histograms
      for(unsigned int ix=1;ix<3;++ix) {
	stringstream ss;
	ss << "TMP/n" << ix;
	book(_nMeson[ix], ss.str());
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      if(ntotal!=3) vetoEvent;
      if(nCount[310]==1 &&
	 ((nCount[ 211]==1&&nCount[-321]==1)||
	  (nCount[-211]==1&&nCount[ 321]==1)))
	 _nMeson[1]->fill();
      else if(nCount[321]==1 &&
	      nCount[-321]==1 && nCount[111]==1)
	 _nMeson[2]->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<3;++ix) {
	double sigma = _nMeson[ix]->val();
	double error = _nMeson[ix]->err();
    	sigma *= crossSection()/ sumOfWeights() /nanobarn;
    	error *= crossSection()/ sumOfWeights() /nanobarn; 
	Scatter2D temphisto(refData(ix, 1, 1));
    	Scatter2DPtr mult;
    	book(mult, ix, 1, 1);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS()/MeV, x-ex2.first, x+ex2.second)) {
	    mult->addPoint(x, sigma, ex, make_pair(error,error));
	  }
	  else {
	    mult->addPoint(x, 0., ex, make_pair(0.,.0));
	  }
	}
      }
    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _nMeson[3];
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(DM2_1991_I318558);


}
