// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BABAR_2006_I731865 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2006_I731865);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      _nPhiEta =  bookCounter("TMP/phieta");
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
      foreach (const Particle& p, ufs.particles()) {
	if(p.children().empty()) continue;
	if(p.pdgId()!=333) continue;
	map<long,int> nRes=nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	foreach (const Particle& p2, ufs.particles()) {
	  if(p2.pdgId()!=221) continue;
	  if(p2.parents()[0].isSame(p)) continue;
	  map<long,int> nResB = nRes;
	  int ncountB = ncount;
	  findChildren(p2,nResB,ncountB);
	  if(ncountB!=0) continue;
	  bool matched2 = true;
	  for(auto const & val : nResB) {
	    if(val.second!=0) {
	      matched2 = false;
	      break;
	    }
	  }
	  if(matched2) {
	    _nPhiEta->fill(event.weight());
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double sigma = _nPhiEta->val();
      double error = _nPhiEta->err();
      sigma *= crossSection()/ sumOfWeights() /femtobarn;
      error *= crossSection()/ sumOfWeights() /femtobarn;
      
      Scatter2D temphisto(refData(1, 1, 2));
      Scatter2DPtr  mult = bookScatter2D(1, 1, 2);
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
    CounterPtr _nPhiEta;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BABAR_2006_I731865);

}
