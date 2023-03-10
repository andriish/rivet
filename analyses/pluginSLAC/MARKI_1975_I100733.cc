// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MARKI_1975_I100733 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKI_1975_I100733);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(FinalState(), "FS");
      // Book histograms
      book(_c_hadrons, "/TMP/sigma_hadrons");
      book(_c_muons, "/TMP/sigma_muons");
      if(isCompatibleWithSqrtS(3*GeV))
      	book(_h_charged, 3, 1, 1);
      else if(isCompatibleWithSqrtS(4.8*GeV))
        book(_h_charged, 3, 1, 2);
      else if(isCompatibleWithSqrtS(7.4*GeV))
        book(_h_charged, 3, 1, 3);
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
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 &&
      	 ntotal==2+nCount[22])
      	_c_muons->fill();
      // everything else
      else {
      	_c_hadrons->fill();
	if(_h_charged) {
	  for (const Particle& p : fs.particles()) {
	    if(PID::isCharged(p.pid())) {
	      double x = 2.*p.p3().mod()/sqrtS();
	      _h_charged->fill(x);
	    }
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      if(_h_charged) {
      	scale(_h_charged, crossSection()/ sumOfWeights() /microbarn*sqr(sqrtS()));
      }
      // R
      Scatter1D R = *_c_hadrons/ *_c_muons;
      double              rval = R.point(0).x();
      pair<double,double> rerr = R.point(0).xErrs();
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      double sig_h = _c_hadrons->val()*fact;
      double err_h = _c_hadrons->err()*fact;
      double sig_m = _c_muons  ->val()*fact;
      double err_m = _c_muons  ->err()*fact;
      Scatter2D temphisto(refData(1, 1, 1));
      Scatter2DPtr hadrons;
      book(hadrons, 1,1,1);
      Scatter2DPtr muons;
      book(muons, "sigma_muons"  );
      Scatter2DPtr mult;
      book(mult, 2,1,1);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
      	const double x  = temphisto.point(b).x();
      	pair<double,double> ex = temphisto.point(b).xErrs();
      	pair<double,double> ex2 = ex;
      	if(ex2.first ==0.) ex2. first=0.0001;
      	if(ex2.second==0.) ex2.second=0.0001;
      	if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
      	  mult   ->addPoint(x, rval, ex, rerr);
      	  hadrons->addPoint(x, sig_h, ex, make_pair(err_h,err_h));
      	  muons  ->addPoint(x, sig_m, ex, make_pair(err_m,err_m));
      	}
      	else {
      	  mult   ->addPoint(x, 0., ex, make_pair(0.,.0));
      	  hadrons->addPoint(x, 0., ex, make_pair(0.,.0));
      	  muons  ->addPoint(x, 0., ex, make_pair(0.,.0));
      	}
      }
    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _c_hadrons, _c_muons;
    Histo1DPtr _h_charged;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MARKI_1975_I100733);


}
