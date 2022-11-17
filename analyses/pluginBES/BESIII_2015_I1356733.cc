// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief 
  class BESIII_2015_I1356733 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2015_I1356733);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<double>("PID", 20445);
      // projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid== 20443 ||
				Cuts::pid==100443), "UFS");
      // histograms
      for(unsigned int ix=0;ix<2;++ix)
	book(_c[ix],"TMP/c_"+toString(ix));
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
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      for(const Particle & psi : apply<UnstableParticles>(event, "UFS").particles()) {
	if(psi.children().empty()) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(psi,nRes,ncount);
	int nphoton = psi.pid()==100443 ? 0 : 1;
	if(ncount!=2+nphoton) continue;
	bool matched = true;
	for(auto const & val : nRes) {
	  if(abs(val.first)==211) {
	    if(val.second !=1) {
	      matched = false;
	      break;
	    }
	  }
	  else if(val.first==22) {
	    if(val.second!=nphoton) {
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
	  if(nphoton==0) {
	    _c[1]->fill();
	    break;
	  }
	  else {
	    Particle parent = psi.parents()[0];
	    if(parent.pid()==_pid && parent.children().size()==2 &&
	       (parent.children()[0].pid()==22 ||
		parent.children()[1].pid()==22 )) {
	      _c[0]->fill();
	      break;
	    }
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /picobarn; 
      for(unsigned int ix=0;ix<2;++ix) {
	double sigma = _c[ix]->val()*fact;
	double error = _c[ix]->err()*fact;
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
    CounterPtr _c[2];
    int _pid;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2015_I1356733);

}
