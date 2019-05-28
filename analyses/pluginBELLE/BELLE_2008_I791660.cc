// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BELLE_2008_I791660 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BELLE_2008_I791660);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      _nLambda = bookCounter("/TMP/nLambda");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      foreach(const Particle &child, p.children()) {
	if(child.children().empty()) {
	  nRes[child.pdgId()]-=1;
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
      foreach (const Particle& p, fs.particles()) {
	nCount[p.pdgId()] += 1;
	++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");

      for(unsigned int ix=0;ix<ufs.particles().size();++ix) {
	const Particle& p1 = ufs.particles()[ix];
	if(abs(p1.pdgId())!=4122) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p1,nRes,ncount);
	bool matched=false;
	for(unsigned int iy=0;iy<ufs.particles().size();++iy) {
	  if(ix==iy) continue;
	  const Particle& p2 = ufs.particles()[iy];
	  if(p2.pdgId() != -p1.pdgId()) continue;
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
	    break;
	  }
	}
	if(matched) {
	  _nLambda->fill(event.weight());
	  break;
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double sigma = _nLambda->val();
      double error = _nLambda->err();
      sigma *= crossSection()/ sumOfWeights() /nanobarn;
      error *= crossSection()/ sumOfWeights() /nanobarn; 
      Scatter2D temphisto(refData(1, 1, 1));
      Scatter2DPtr  mult = bookScatter2D(1, 1, 1);
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

    //@}


    /// @name Histograms
    //@{
    CounterPtr _nLambda;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BELLE_2008_I791660);


}
