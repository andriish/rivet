// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class SND_2001_I533574 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(SND_2001_I533574);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      _nKpKm = bookCounter("TMP/KpKm");
      _nK0K0 = bookCounter("TMP/K0K0");
      _n3pi  = bookCounter("TMP/3pi");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0);
      foreach (const Particle& p, fs.particles()) {
	nCount[p.pdgId()] += 1;
	++ntotal;
      }
      if(ntotal==2) {
	if(nCount[321]==1 && nCount[-321]==1)
	  _nKpKm->fill(event.weight());
	else if(nCount[130]==1 && nCount[310]==1)
	  _nK0K0->fill(event.weight());
      }
      else if(ntotal==3 && nCount[211] == 1 && nCount[-211] == 1 && nCount[111] == 1)
	_n3pi->fill(event.weight());

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int iy=1;iy<3;++iy) {
	for(unsigned int ix=1;ix<5;++ix) {
	  double sigma,error;
	  if(ix==1) {
	    sigma = _nKpKm->val();
	    error = _nKpKm->err();
	  }
	  else if(ix==2) {
	    sigma = _nK0K0->val();
	    error = _nK0K0->err();
	  }
	  else if(ix==3) {
	    sigma = _nK0K0->val();
	    error = _nK0K0->err();
	  }
	  else if(ix==4) {
	    sigma = _n3pi->val();
	    error = _n3pi->err();
	  }
	  sigma *= crossSection()/ sumOfWeights() /nanobarn;
	  error *= crossSection()/ sumOfWeights() /nanobarn; 
	  Scatter2D temphisto(refData(iy, 1, ix));
	  Scatter2DPtr  mult = bookScatter2D(iy, 1, ix);
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
    }
    
    //@}


    /// @name Histograms
    //@{
    CounterPtr _nKpKm,_nK0K0,_n3pi;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(SND_2001_I533574);


}
