// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BABAR_2005_I676691 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2005_I676691);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      book(_c2pip2pim  , "TMP/2pip2pim");
      book(_cKpKmpippim, "TMP/KpKmpippim");
      book(_c2Kp2Km    , "TMP/2Kp2Km");

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

      if(ntotal!=4) vetoEvent;

      if( nCount[211]==2 && nCount[-211]==2)
	_c2pip2pim->fill();
      else if(nCount[321]==1 && nCount[-321]==1 && nCount[211]==1 && nCount[-211]==1)
	_cKpKmpippim->fill();
      else if( nCount[321]==2 && nCount[-321]==2)
	_c2Kp2Km->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      for(unsigned int ix=1;ix<4;++ix) {
	double sigma = 0., error = 0.;
	if(ix==1) {
	  sigma =  _c2pip2pim->val();
	  error =  _c2pip2pim->err();
	}
	else if(ix==2) {
	  sigma = _cKpKmpippim->val();
	  error = _cKpKmpippim->err();
	}
     	else if(ix==3) {
	  sigma =  _c2Kp2Km->val();
	  error =  _c2Kp2Km->err();
	}
    	sigma *= crossSection()/ sumOfWeights() /nanobarn;
    	error *= crossSection()/ sumOfWeights() /nanobarn;
	Scatter2D temphisto(refData(ix, 1, 1));
	Scatter2DPtr  mult;
        book(mult, ix, 1, 1);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
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
    CounterPtr _c2pip2pim, _cKpKmpippim, _c2Kp2Km;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BABAR_2005_I676691);


}
