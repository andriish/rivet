// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e-> KS0 Kpi pi+pi-
  class CMD3_2022_I2108984 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD3_2022_I2108984);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==20223), "UFS");
      // counters
      for(unsigned int ix=0;ix<3;++ix)
	book(_c[ix]  , "/TMP/c"+toString(ix+1));
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for(const Particle &child : p.children()) {
	if(child.children().empty()) {
	  --nRes[child.pid()];
	  --ncount;
	}
	else
	  findChildren(child,nRes,ncount);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find the final-state particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // the FS cross section
      bool foundFS = false;
      if(ntotal==5 && nCount[310]==1 &&
	 ( ( nCount[211]==2 && nCount[-211]==1 && nCount[-321] ==1 ) ||
	   ( nCount[211]==1 && nCount[-211]==2 && nCount[ 321] ==1 ) ) ) {
	foundFS = true;
	_c[0]->fill();
      }
      // check for intermediate f_1
      for (const Particle& p : apply<FinalState>(event, "UFS").particles()) {
	if(p.children().empty()) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	if(ncount!=2) continue;
	bool matched = true;
	for(auto const & val : nRes) {
	  if(abs(val.first)==211 ) {
	    if(val.second !=1) {
	      matched = false;
	      break;
	    }
	  }
	  else if(val.second!=0) {
	    matched = false;
	    break;
	  }
	}
	if (matched) {
	  _c[2]->fill();
	  if (foundFS) _c[1]->fill();
	  break;
	}
      }


      
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ix=1;ix<3;++ix) {
	for(unsigned int iy=1;iy<ix+1;++iy) {
	  double sigma = _c[ix+iy-2]->val()*fact;
	  double error = _c[ix+iy-2]->err()*fact;
	  Scatter2D temphisto(refData(ix, 1, iy));
	  Scatter2DPtr  mult;
	  book(mult, ix, 1, iy);
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

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMD3_2022_I2108984);

}
