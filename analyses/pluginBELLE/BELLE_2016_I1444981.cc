// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief gamma gamma -> p K+ pbar K-
  class BELLE_2016_I1444981 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2016_I1444981);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::abspid==101234), "UFS");
      // counters
      for(unsigned int ix=0;ix<2;++ix)
	book(_c[ix],"TMP/c_"+toString(ix+1));
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
      // find the final-state particles
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // p pbar K+ K-
      if(ntotal==4 &&
	 nCount[2212]==1 && nCount[-2212]==1 &&
	 nCount[ 321]==1 && nCount[ -321]==1) _c[0]->fill();
      // intermediate Lambda(1520)
      for (const Particle & lam : apply<UnstableParticles>(event, "UFS").particles()) {
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(lam,nRes,ncount);
	int sign = lam.pid()/lam.abspid();
	if(ncount!=2) continue;
	bool matched = true;
	for(auto const & val : nRes) {
	  if(val.first==-sign*2212 ||
	     val.first==-sign*321) {
	    if(val.second!=0) {
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
	  _c[1]->fill();
	  break;
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/picobarn/sumOfWeights();
      // loop over tables in paper
      for(unsigned int ix=0;ix<2;++ix) {
	double sigma = _c[ix]->val()*fact;
	double error = _c[ix]->err()*fact;
	Scatter2D temphisto(refData(ix+1, 1, 1));
	Scatter2DPtr mult;
	book(mult, ix+1, 1, 1);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS(), x-ex2.first, x+ex2.second)) {
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
    CounterPtr _c[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2016_I1444981);

}
