// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMD2_2000_I532970 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMD2_2000_I532970);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      _numOmegaPiPi = bookCounter("TMP/OmegaPiPi");
      _numEtaPiPi   = bookCounter("TMP/EtaPiPi");
      _num5Pi       = bookCounter("TMP/5Pi");

    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      foreach(const Particle &child, p.children()) {
	if(child.children().empty()) {
	  --nRes[child.pdgId()];
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
      bool foundRes = false;
      foreach (const Particle& p, ufs.particles()) {
	if(p.children().empty()) continue;
	// find the eta
	if(p.pdgId()==221) {
	  map<long,int> nRes = nCount;
	  int ncount = ntotal;
	  findChildren(p,nRes,ncount);
	  // eta pi+pi-
	  if(ncount!=2) continue;
	  bool matched = true;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==211) {
	      if(val.second !=1) {
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
	    _numEtaPiPi->fill(event.weight());
	    foundRes = true;
	  }
	}
	// find the eta
	else if(p.pdgId()==223) {
	  map<long,int> nRes = nCount;
	  int ncount = ntotal;
	  findChildren(p,nRes,ncount);
	  // eta pi+pi-
	  if(ncount!=2) continue;
	  bool matched = true;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==211) {
	      if(val.second !=1) {
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
	    _numOmegaPiPi->fill(event.weight());
	    foundRes = true;
	  }
	}
      }

      if(foundRes) vetoEvent;
      if(nCount[-211]==2&&nCount[211]==2&&nCount[111]==1)
	_num5Pi->fill(event.weight());
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<4;++ix) {
	double sigma,error;
	if(ix==1) {
	  sigma = _numOmegaPiPi->val();
	  error = _numOmegaPiPi->err();
	}
	else if(ix==2) {
	  sigma = _numEtaPiPi->val();
	  error = _numEtaPiPi->err();
	}
	else {
	  sigma = _num5Pi->val();
	  error = _num5Pi->err();
	}
	sigma *= crossSection()/ sumOfWeights() /nanobarn;
	error *= crossSection()/ sumOfWeights() /nanobarn; 
	Scatter2D temphisto(refData(ix, 1, 1));
	Scatter2DPtr  mult = bookScatter2D(ix, 1, 1);
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

    //@}

    /// @name Histograms
    //@{
    CounterPtr _numEtaPiPi,_numOmegaPiPi,_num5Pi;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMD2_2000_I532970);


}
