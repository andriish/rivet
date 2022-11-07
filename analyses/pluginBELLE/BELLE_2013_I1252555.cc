// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > omega pi0, K* Kbar K2 Kbar
  class BELLE_2013_I1252555 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2013_I1252555);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // counters
      for(unsigned int ix=0;ix<5;++ix)
	book(_c[ix],"TMP/c_"+toString(ix+1));
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for(const Particle &child : p.children()) {
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
      // final state particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // resonances
      for(const Particle & p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==223 or
										Cuts::abspid==313 or
										Cuts::abspid==323 or
										Cuts::abspid==315 or
										Cuts::abspid==325)) {
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	if(ncount!=1) continue;
	unsigned int imode=0;
	long pid = 111;
	if(p.abspid()==313) {
	  imode=1;
	  pid=310;
	}
	else if(p.abspid()==323) {
	  imode=2;
	  pid = p.pid()>0 ? -321 : 321;
	}
	else if(p.abspid()==315) {
	  imode=3;
	  pid=310;
	}
	else if(p.abspid()==325) {
	  imode=4;
	  pid = p.pid()>0 ? -321 : 321;
	}
	bool matched=true;
	for(auto const & val : nRes) {
	  if((pid==310 && (val.first==310 || val.second==130)) ||
	     val.first==pid) {
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
	  _c[imode]->fill();
	  break;
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /femtobarn;
      for(unsigned int iy=0;iy<5;++iy) {
	double sigma = _c[iy]->val()*fact;
	double error = _c[iy]->err()*fact;
	Scatter2D temphisto(refData( 1+iy, 1, 1));
    	Scatter2DPtr  mult;
        book(mult, 1+iy, 1, 1);
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
    CounterPtr _c[5];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2013_I1252555);

}
