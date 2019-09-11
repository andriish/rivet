// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class BABAR_2017_I1621593 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2017_I1621593);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      _num4pi   = bookCounter("TMP/num4");
      _numOmega = bookCounter("TMP/numOmega");

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
      if(ntotal!=4) vetoEvent;
      if(nCount[-211]==1&&nCount[211]==1&&nCount[111]==2) {
	_num4pi->fill(event.weight());
	const FinalState& ufs = apply<FinalState>(event, "UFS");
	if(!ufs.particles(Cuts::pid==223).empty()) 
	  _numOmega->fill(event.weight());
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      for(unsigned int ix=1;ix<3;++ix) {
	double sigma,error;
	if(ix==1) {
	  sigma = _num4pi->val();
	  error = _num4pi->err();
	}
	else {
	  sigma = _numOmega->val();
	  error = _numOmega->err();
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
    CounterPtr _num4pi,_numOmega;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BABAR_2017_I1621593);


}
