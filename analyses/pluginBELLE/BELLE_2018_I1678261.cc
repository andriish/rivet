// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief 
  class BELLE_2018_I1678261 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2018_I1678261);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==20553 or
				Cuts::pid==555), "UFS");
      book(_c,"TMP/nChi");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
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
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p: fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      // loop over any chi mesons
      for(const Particle & chi : ufs.particles()) {
	bool matched = false;
	if(chi.children().empty()) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(chi,nRes,ncount);
	matched = true;
	for(auto const & val : nRes) {
	  if(abs(val.first)==PID::PIPLUS || val.first==PID::PI0) {
	    if (val.second!=1) {
	      matched = false;
	      break;
	    }
	  }
	  else if (val.second!=0) {
	    matched = false;
	    break;
	  }
	}
	if(!matched) continue;
	_c->fill();
	break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact =  crossSection()/ sumOfWeights() /picobarn;
      double sigma = _c->val()*fact;
      double error = _c->err()*fact;
      Scatter2D temphisto(refData(1, 1, 1));
      Scatter2DPtr  mult;
      book(mult, 1, 1, 1);
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

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2018_I1678261);

}