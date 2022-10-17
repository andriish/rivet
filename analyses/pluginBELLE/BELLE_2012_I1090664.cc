// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief 
  class BELLE_2012_I1090664 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2012_I1090664);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==223 or Cuts::pid==333), "UFS");
      // counters
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
      // find the final-state particles
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // find any rho mesons
      Particles mesons=apply<UnstableParticles>(event, "UFS").particles();
      for (unsigned int ix=0;ix<mesons.size();++ix) {
       	if(mesons[ix].children().empty()) continue;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(mesons[ix],nRes,ncount);
      	bool matched = false;
      	for (unsigned int iy=ix+1;iy<mesons.size();++iy) {
      	  if(mesons[iy].children().empty()) continue;
       	  map<long,int> nRes2=nRes;
       	  int ncount2 = ncount;
       	  findChildren(mesons[iy],nRes2,ncount2);
       	  if(ncount2 !=0 ) continue;
       	  matched=true;
      	  for(auto const & val : nRes2) {
      	    if(val.second!=0) {
      	      matched = false;
      	      break;
      	    }
      	  }
      	  if(matched) {
	  if(mesons[ix].pid()!=mesons[iy].pid())
	    _c[0]->fill();
	  else if(mesons[ix].pid()==333)
	    _c[1]->fill();
	  else if(mesons[ix].pid()==223)
	    _c[2]->fill();
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
      for(unsigned int ix=0;ix<3;++ix) {
	double sigma = _c[ix]->val()*fact;
	double error = _c[ix]->err()*fact;
	for(unsigned int iy=0;iy<2;++iy) {
	  Scatter2D temphisto(refData(1, ix+1, 1+iy));
	  Scatter2DPtr mult;
	  book(mult, 1, ix+1, 1+iy);
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
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2012_I1090664);

}
