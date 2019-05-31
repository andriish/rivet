// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class BABAR_2014_I1287920 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2014_I1287920);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      _nKSKL     = bookCounter("TMP/nKSKL");
      _nKSKLpipi = bookCounter("TMP/nKSKLpipi");
      _nKSKSpipi = bookCounter("TMP/nKSKSpipi");
      _nKSKSKpKm = bookCounter("TMP/nKSKSKpKm");
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
      
      if(ntotal==2 && nCount[130]==1 && nCount[310]==1)
	_nKSKL->fill(event.weight());
      else if( ntotal==4 &&  nCount[130]==1 && nCount[310]==1 && nCount[211]==1 && nCount[-211]==1)
	_nKSKLpipi->fill(event.weight());
      else if( ntotal==4 && nCount[310]==2 && nCount[211]==1 && nCount[-211]==1 )
	_nKSKSpipi->fill(event.weight());
      else if( ntotal==4 && nCount[310]==2 && nCount[321]==1 && nCount[-321]==1)
	_nKSKSKpKm->fill(event.weight());

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      for(unsigned int ix=9;ix<13;++ix) {
	double sigma,error;
	if(ix==9) {
	  sigma =  _nKSKL->val();
	  error =  _nKSKL->err();
	}
	else if(ix==10) {
	  sigma = _nKSKLpipi->val();
	  error = _nKSKLpipi->err();
	}
     	else if(ix==11) {
	  sigma = _nKSKSpipi->val();
	  error = _nKSKSpipi->err();
	}
     	else if(ix==12) {
	  sigma = _nKSKSKpKm->val();
	  error = _nKSKSKpKm->err();
	}
    	sigma *= crossSection()/ sumOfWeights() /nanobarn;
    	error *= crossSection()/ sumOfWeights() /nanobarn;
	Scatter2D temphisto(refData(ix, 1, 1));
	Scatter2DPtr  mult = bookScatter2D(ix, 1, 1);
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
    CounterPtr _nKSKL,_nKSKLpipi,_nKSKSpipi,_nKSKSKpKm;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BABAR_2014_I1287920);


}
