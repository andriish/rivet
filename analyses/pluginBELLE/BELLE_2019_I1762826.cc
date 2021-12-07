// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  e+e- -> Ds D_s1
  class BELLE_2019_I1762826 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2019_I1762826);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_c_Ds    ,"/TMP/c_Ds"    );
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
      const FinalState& fs = apply<FinalState>(event, "FS");
      // total analyse final state
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      // unstable charm analysis
      Particles ds = apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==431 or Cuts::abspid==10433);
      for(unsigned int ix=0;ix<ds.size();++ix) {
	const Particle& p1 = ds[ix];
       	int id1 = abs(p1.pid());
	// check fs
	bool fs = true;
	for (const Particle & child : p1.children()) {
	  if(child.pid()==p1.pid()) {
	    fs = false;
	    break;
	  }
	}
      	if(!fs) continue;
      	// find the children
      	map<long,int> nRes = nCount;
      	int ncount = ntotal;
      	findChildren(p1,nRes,ncount);
      	bool matched=false;
       	int sign = p1.pid()/id1;
       	for(unsigned int iy=ix+1;iy<ds.size();++iy) {
      	  const Particle& p2 = ds[iy];
      	  fs = true;
      	  for (const Particle & child : p2.children()) {
      	    if(child.pid()==p2.pid()) {
      	      fs = false;
      	      break;
      	    }
      	  }
      	  if(!fs) continue;
      	  if(p2.pid()/abs(p2.pid())==sign) continue;
	  int id2 = abs(p2.pid());
	  if(!p2.parents().empty() && p2.parents()[0].pid()==p1.pid())
	    continue;
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
	  if(matched) {
       	    if((id1==431 && id2==10433) || (id1==10433 && id2==431)) {
	      Particle Ds2 = id1==10433 ? p1 : p2;
	      matched = false;
	      if(Ds2.children().size()!=2) continue;
	      if(Ds2.pid()==10433 &&
		 ((Ds2.children()[0].pid()==423 && Ds2.children()[1].pid()==321) ||
		  (Ds2.children()[1].pid()==423 && Ds2.children()[0].pid()==321) ))
		matched = true;
	      else if (Ds2.pid()==-10433 &&
		 ((Ds2.children()[0].pid()==-423 && Ds2.children()[1].pid()==-321) ||
		  (Ds2.children()[1].pid()==-423 && Ds2.children()[0].pid()==-321) ))
		matched = true;
	      if(matched) {
		_c_Ds->fill();
		break;
	      }
	    }
	  }
	}
	if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights()/picobarn;
      double sigma = _c_Ds->val()*fact;
      double error = _c_Ds->err()*fact;
      Scatter2D temphisto(refData(1, 1, 1));
      Scatter2DPtr     mult;
      book(mult, 1, 1, 1);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
	const double x  = temphisto.point(b).x();
	pair<double,double> ex = temphisto.point(b).xErrs();
	pair<double,double> ex2 = ex;
	if(ex2.first ==0.) ex2. first=0.0001;
	if(ex2.second==0.) ex2.second=0.0001;
	if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
	  mult   ->addPoint(x, sigma, ex, make_pair(error,error));
	}
	else {
	  mult   ->addPoint(x, 0., ex, make_pair(0.,.0));
	}
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _c_Ds;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2019_I1762826);

}
