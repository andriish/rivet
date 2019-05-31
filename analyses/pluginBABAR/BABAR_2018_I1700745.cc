// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BABAR_2018_I1700745 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2018_I1700745);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {


      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      _num5Pi        = bookCounter("TMP/5Pi");
      _num2PiEta     = bookCounter("TMP/2PiEta");
      _numOmegaPiPi  = bookCounter("TMP/OmegaPiPi");
      _num4PiEta     = bookCounter("TMP/4PiEta");
      _numOmegaPiEta = bookCounter("TMP/OmegaPiEta");
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
      if(ntotal==5 && nCount[211]==1 && nCount[-211]==1 && nCount[111]==3)
	_num5Pi->fill(event.weight());
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      foreach (const Particle& p, ufs.particles()) {
     	if(p.children().empty()) continue;
	// find eta/omegas
      	if(p.pdgId()==221 || p.pdgId()==223 ) {
       	  map<long,int> nRes = nCount;
       	  int ncount = ntotal;
       	  findChildren(p,nRes,ncount);
	  // eta
	  if(p.pdgId()==221) {
	    // 2 pi eta
	    if(ncount==2) {
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
	      if(matched)
		_num2PiEta->fill(event.weight());
	    }
	    // 4 pi eta
	    else if(ncount==4) {
	      bool matched = true;
	      for(auto const & val : nRes) {
		if(abs(val.first)==211) {
		  if(val.second !=1) {
		    matched = false;
		    break;
		  }
		}
		else if(abs(val.first)==111) {
		  if(val.second !=2) {
		    matched = false;
		    break;
		  }
		}
		else if(val.second!=0) {
		  matched = false;
		  break;
		}
	      }
	      if(matched)
		_num4PiEta->fill(event.weight());
	    }
	    // pi0 omega eta
	    foreach (const Particle& p2, ufs.particles()) {
	      if(p2.pdgId()!=223) continue;
	      map<long,int> nResB = nRes;
	      int ncountB = ncount;
	      findChildren(p2,nResB,ncountB);
	      if(ncountB!=1) continue;
	      bool matched = true;
	      for(auto const & val : nResB) {
		if(abs(val.first)==111) {
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
		_numOmegaPiEta->fill(event.weight());
		break;
	      }
	    }
	  }
	  else {
	    if(ncount!=2) continue;	    
	    bool matched = true;
	    for(auto const & val : nRes) {
	      if(abs(val.first)==111) {
		if(val.second !=2) {
		  matched = false;
		  break;
		}
	      }
	      else if(val.second!=0) {
		matched = false;
		break;
	      }
	    }
	    if(matched)
	      _numOmegaPiPi->fill(event.weight());
	  }
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<6;++ix) {
	double sigma,error;
	if(ix==1) {
	  sigma = _num5Pi->val();
	  error = _num5Pi->err();
	}
	else if(ix==2) {
	  sigma = _num2PiEta->val();
	  error = _num2PiEta->err();
	}
	else if(ix==3) {
	  sigma = _numOmegaPiPi->val();
	  error = _numOmegaPiPi->err();
	}
	else if(ix==4) {
	  sigma = _num4PiEta->val();
	  error = _num4PiEta->err();
	}
	else if(ix==5) {
	  sigma = _numOmegaPiEta->val();
	  error = _numOmegaPiEta->err();
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
    CounterPtr _num5Pi,_num2PiEta,_numOmegaPiPi,
      _num4PiEta,_numOmegaPiEta;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BABAR_2018_I1700745);


}
