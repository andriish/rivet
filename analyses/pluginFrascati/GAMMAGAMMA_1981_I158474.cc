// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class GAMMAGAMMA_1981_I158474 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(GAMMAGAMMA_1981_I158474);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      _n3pi  = bookCounter("TMP/n3pi");
      _n4pi  = bookCounter("TMP/n4pi");
      _n5pi  = bookCounter("TMP/n5pi");
      _n6pi  = bookCounter("TMP/n6pi");
      _n35pi = bookCounter("TMP/n35pi");
      _n46pi = bookCounter("TMP/n46pi");
      _nC2   = bookCounter("TMP/nC2");
      _nC4   = bookCounter("TMP/nC4");
      _nmu   = bookCounter("TMP/nmu");
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
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 &&
	 ntotal==2+nCount[22])
	_nmu->fill(event.weight());
      else {
	if(ntotal==3 && nCount[211] == 1 && nCount[-211]==1 && nCount[111]==1 ) {
	  _n3pi->fill(event.weight());
	}
	if(ntotal==4 && nCount[211] == 1 && nCount[-211]==1 && nCount[111]==2 ) {
	  _n4pi->fill(event.weight());
	}
	if(ntotal==5 && nCount[211] == 2 && nCount[-211]==2 && nCount[111]==1 ) {
	  _n5pi->fill(event.weight());
	}
	if(ntotal==6 && nCount[211] == 2 && nCount[-211]==2 && nCount[111]==2 ) {
	  _n6pi->fill(event.weight());
	}
	if(nCount[211] == 1 && nCount[-211]==1 && ntotal == 2+nCount[111]) {
	  _nC2->fill(event.weight());
	}
	if(nCount[211] == 2 && nCount[-211]==2 && ntotal == 4+nCount[111]) {
	  _nC4->fill(event.weight());
	}
	if((nCount[211]+nCount[-211]+nCount[111])==ntotal ) {
	  if(ntotal==3 || ntotal ==5)
	    _n35pi->fill(event.weight());
	  else if(ntotal==4 || ntotal==6) 
	    _n46pi ->fill(event.weight());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      for(unsigned int ix=1;ix<7;++ix) {
	double sigma,error;
	if(ix==1) {
	  sigma = _n3pi->val()*fact;
	  error = _n3pi->err()*fact;
	}
	else if(ix==2) {
	  sigma = _n4pi->val()*fact;
	  error = _n4pi->err()*fact;
	}
	else if(ix==3) {
	  sigma = _n5pi->val()*fact;
	  error = _n5pi->err()*fact;
	}
	else if(ix==4) {
	  sigma = _n6pi->val()*fact;
	  error = _n6pi->err()*fact;
	}
	else if(ix==5) {
	  sigma = _n35pi->val()*fact;
	  error = _n35pi->err()*fact;
	}
	else if(ix==6) {
	  sigma = _n46pi->val()*fact;
	  error = _n46pi->err()*fact;
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
	    mult->addPoint(x, sigma, ex, make_pair(error,error));
	  }
	  else {
	  mult->addPoint(x, 0., ex, make_pair(0.,.0));
	  }
	}
      }
      for(unsigned int ix=1;ix<3;++ix) {
	Scatter1D R = (ix==1? *_nC2 : *_nC4)/ *_nmu;
	double              rval = R.point(0).x();
	pair<double,double> rerr = R.point(0).xErrs();
	double sig_h = (ix ==1 ? _nC2 : _nC4)->val()*fact;
	double err_h = (ix ==1 ? _nC2 : _nC4)->err()*fact;
	double sig_m = _nmu->val()*fact;
	double err_m = _nmu->err()*fact;
	Scatter2D temphisto(refData(2, 1, ix));
	ostringstream title;
	if(ix==1)
	  title << "sigma_2pi";
	else
	  title << "sigma_4pi";
	Scatter2DPtr hadrons  = bookScatter2D(title.str());
	Scatter2DPtr muons;
	if(ix==1) muons = bookScatter2D("sigma_muons");
	Scatter2DPtr     mult = bookScatter2D(2,1,ix);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
	    mult   ->addPoint(x, rval, ex, rerr);
	    hadrons->addPoint(x, sig_h, ex, make_pair(err_h,err_h));
	    if(ix==1) muons  ->addPoint(x, sig_m, ex, make_pair(err_m,err_m));
	  }
	  else {
	    mult   ->addPoint(x, 0., ex, make_pair(0.,.0));
	    hadrons->addPoint(x, 0., ex, make_pair(0.,.0));
	    if(ix==1) muons  ->addPoint(x, 0., ex, make_pair(0.,.0));
	  }
	}
      }
    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _n3pi,_n4pi,_n5pi,_n6pi,_n35pi,_n46pi,_nC2,_nC4,_nmu;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(GAMMAGAMMA_1981_I158474);


}
