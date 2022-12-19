// -*- C++ -*-
#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief e+ e- > D Dbar
  class BESII_2008_I802372 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2008_I802372);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::abspid==411 ||
				Cuts::abspid==421), "UFS");
      // histos
      for(unsigned int ix=0;ix<3;++ix)
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

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for(const Particle & p1 : ufs.particles()) {
	if(p1.pid()<0) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p1,nRes,ncount);
	bool matched=false;
	for(const Particle & p2 : ufs.particles()) {
     	  if(p2.pid()!=-p1.pid()) continue;
   	  map<long,int> nRes2 = nRes;
   	  int ncount2 = ncount;
   	  findChildren(p2,nRes2,ncount2);
   	  if(ncount2!=0) continue;
   	  matched=true;
   	  for(auto const & val : nRes2) {
   	    if(val.second!=0) {
   	      matched = false;
   	      break;
   	    }
   	  }
   	  if(matched) break;
   	}
   	if(matched) {
	  _c[2]->fill();
	  if(p1.abspid()==421)      _c[0]->fill();
	  else if(p1.abspid()==411) _c[1]->fill();
	  break;
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix) {
       	double sigma = _c[ix]->val() * crossSection()/ sumOfWeights() /nanobarn;
	double error = _c[ix]->err() * crossSection()/ sumOfWeights() /nanobarn;
       	Scatter2D temphisto(refData(1, 1, 1+ix));
       	Scatter2DPtr  mult;
	book(mult, 1, 1, 1+ix);
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
    CounterPtr _c[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2008_I802372);

}
