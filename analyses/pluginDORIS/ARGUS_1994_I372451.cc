// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief gamma gamma -> phi rho or omega
  class ARGUS_1994_I372451 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1994_I372451);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // book histos
      if(inRange(sqrtS()/GeV,1.5,3.5)) {
	for(unsigned int ix=0;ix<2;++ix)
	  book(_nMeson[ix],"TMP/nMeson_"+toString(ix+1));
      }
      else
	throw Error("Invalid CMS energy for ARGUS_1994_I372451");
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
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles phis = ufs.particles(Cuts::pid==333);
      Particles rhos = ufs.particles(Cuts::pid==113 or Cuts::pid==223);
      for(const Particle & phi : phis) {
	bool matched=false;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(phi,nRes,ncount);
	for(const Particle & rho : rhos) {
	  map<long,int> nRes2=nRes;
	  int ncount2 = ncount;
	  findChildren(rho,nRes2,ncount2);
	  if(ncount2 !=0 ) continue;
	  matched = true;
	  for(auto const & val : nRes2) {
	    if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) {
	    if(rho.pid()==113)
	      _nMeson[0]->fill();
	    else
	      _nMeson[1]->fill();
	    break;
	  }
	}
	if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      // loop over tables in paper
      for(unsigned int ix=0;ix<2;++ix) {
	double sigma = _nMeson[ix]->val()*fact;
	double error = _nMeson[ix]->err()*fact;
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
    CounterPtr _nMeson[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1994_I372451);

}
