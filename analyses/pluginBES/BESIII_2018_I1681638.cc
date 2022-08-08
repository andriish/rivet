// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> p KS0 nbar K-
  class BESIII_2018_I1681638 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1681638);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==101234), "UFS");
      // counters
      for(unsigned int ix=0;ix<3;++ix)
	book(_n[ix]  , "/TMP/c"+toString(ix+1));
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
      // FS
      if(ntotal==4&&nCount[310]==1) {
	if((nCount[ 2212]==1 && nCount[-2112]==1 && nCount[-321]==1) ||
	   (nCount[-2212]==1 && nCount[ 2112]==1 && nCount[ 321]==1))
	  _n[0]->fill();
      }
      for (const Particle& p : apply<FinalState>(event, "UFS").particles()) {
	if(p.children().empty()) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	if(ncount!=2) continue;
	bool matched = true;
	int sign = p.pid()/p.abspid();
	for(auto const & val : nRes) {
	  if(val.first==-sign*2112 ||
	     val.first==310) {
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
	if(matched) {
	  _n[1]->fill();
	  break;
	}
	for(auto const & val : nRes) {
	  if(val.first==-sign*2212 ||
	     val.first==sign*321) {
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
	if(matched) {
	  _n[2]->fill();
	  break;
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/picobarn/sumOfWeights();
      for(unsigned int ix=0;ix<3;++ix) {
	double sigma = _n[ix]->val()*fact;
	double error = _n[ix]->err()*fact;
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
    CounterPtr _n[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1681638);

}
