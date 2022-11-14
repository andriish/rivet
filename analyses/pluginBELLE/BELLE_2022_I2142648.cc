// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  e+e- -> omega chi_b(0,1,2)
  class BELLE_2022_I2142648 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2022_I2142648);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // histos
      for(unsigned int ix=0;ix<3;++ix)
	book(_c[ix],"TMP/nChi_"+toString(ix));
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
      for(const Particle & chi : ufs.particles(Cuts::pid==10551 or
					       Cuts::pid==20553 or
					       Cuts::pid==555)) {
	bool matched = false;
	if(chi.children().empty()) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(chi,nRes,ncount);
	// loop over omega mesons
	for(const Particle & omega : ufs.particles(Cuts::pid==223)) {
	  map<long,int> nRes2 = nRes;
	  int ncount2 = ncount;
	  findChildren(omega,nRes2,ncount2);
	  matched = true;
	  for(auto const & val : nRes2) {
	    if (val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(!matched) continue;
	  if (chi.pid()==10441)     _c[0]->fill();
	  else if(chi.pid()==20443) _c[1]->fill();
	  else if(chi.pid()==445)   _c[2]->fill();
	  break;
	}
	if (matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact =  crossSection()/ sumOfWeights() /picobarn;
      for(unsigned int ix=0;ix<3;++ix) {
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
    CounterPtr _c[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2022_I2142648);

}
