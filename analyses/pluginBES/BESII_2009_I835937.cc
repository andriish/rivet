// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- > light hadrons including KS0
  class BESII_2009_I835937 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2009_I835937);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // histograms
      if (isCompatibleWithSqrtS(3.773))
	book(_h,1,1,2);
      for(unsigned int ix=0;ix<7;++ix) {
	if(ix>0&&ix<5) continue;
	book(_c[ix],"TMP/c_"+toString(ix+1));
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
	if(child.children().empty()) {
	  nRes[child.pid()]-=1;
	  --ncount;
	}
	else
	  findChildren(child,nRes,ncount);
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
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      // first the purely FS particles
      if(ntotal==5 && nCount[310]==1 && ((nCount[-321]==2 && nCount[ 321]==1 && nCount[ 211]==1) ||
					 (nCount[ 321]==2 && nCount[-321]==1 && nCount[-211]==1)))
	_c[0]->fill();
      // loop over unstable particles
      for (const Particle& p : ufs.particles(Cuts::pid==221 || Cuts::pid==113 || Cuts::abspid==213)) {
      	if(p.children().empty()) continue;
      	map<long,int> nRes = nCount;
      	int ncount = ntotal;
      	findChildren(p,nRes,ncount);
	// eta /rho0
	bool matched = false;
	if(p.pid()==221 || p.pid()==113) {
       	  if(ncount==3) {
      	    matched = true;
      	    for(auto const & val : nRes) {
      	      if(abs(val.first)==310) {
      		if(val.second!=1) {
      		  matched = false;
      		  break;
      		}
      	      }
	      else if(abs(val.first)==211) {
      		if(val.second!=1) {
      		  matched = false;
      		  break;
      		}
      	      }
	      else if(abs(val.first)==321) {
      		if(val.second!=1) {
      		  matched = false;
      		  break;
      		}
      	      }
      	      else if(val.second!=0) {
      		matched = false;
      		break;
      	      }

	    }
      	    if(matched) {
	      if(p.pid()==113) _c[6]->fill();
	      else if (p.pid()==221 &&_h) _h->fill(sqrtS());
	    }
	  }
	}
	// rho+/-
	else if(p.abspid()==213) {
	  if(ncount==2) {
	    int sign = p.pid()>0 ? -1 : 1;
	    matched=true;
	    for(auto const & val : nRes) {
      	      if(abs(val.first)==310) {
      		if(val.second!=1) {
      		  matched = false;
      		  break;
      		}
      	      }
	      else if(val.first==321*sign) {
      		if(val.second!=1) {
      		  matched = false;
      		  break;
      		}
      	      }
      	      else if(val.second!=0) {
      		matched = false;
      		break;
      	      }
      	    }
	    if(matched) _c[5]->fill();
	  }
	}
	if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /picobarn;
      if(_h) scale(_h,fact);
      for(unsigned int ix=0;ix<7;++ix) {
	if(ix>0&&ix<5) continue;
	double sigma = _c[ix]->val()*fact;
	double error = _c[ix]->err()*fact;
	Scatter2D temphisto(refData(1, 1, ix+1));
    	Scatter2DPtr  mult;
        book(mult, 1, 1, ix+1);
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

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    CounterPtr _c[7];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2009_I835937);

}
