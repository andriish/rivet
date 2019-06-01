// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ND_1991_I321108 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ND_1991_I321108);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      _nOmegaPi = bookCounter("TMP/OmegaPi");
      _n2Pi     = bookCounter("TMP/2Pi");
      _n3Pi     = bookCounter("TMP/3Pi");
      _n4PiC    = bookCounter("TMP/4PiC");
      _n4PiN    = bookCounter("TMP/4PiN");
      _nEtaPiPi = bookCounter("TMP/EtaPiPi");
      _nKC      = bookCounter("TMP/KC");
      _nKN      = bookCounter("TMP/KN");
      _n5Pi     = bookCounter("TMP/5Pi");
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
      if(ntotal==2) {
	if(nCount[211]==1&&nCount[-211]==1)
	  _n2Pi->fill(event.weight());
	if(nCount[321]==1&&nCount[-321]==1)
	  _nKC->fill(event.weight());
	if(nCount[310]==1&&nCount[130]==1)
	  _nKN->fill(event.weight());
      }
      else if(ntotal==3) {
	if(nCount[211]==1&&nCount[-211]==1&&nCount[111]==1)
	  _n3Pi->fill(event.weight());
      }
      else if(ntotal==4) {
	if(nCount[211]==2&&nCount[-211]==2)
	  _n4PiC->fill(event.weight());
	else if(nCount[211]==1&&nCount[-211]==1&&nCount[111]==2)
	  _n4PiN->fill(event.weight());
      }
      else if(ntotal==5) {
	if(nCount[211]==2&&nCount[-211]==2&&nCount[111]==1)
	  _n5Pi->fill(event.weight());
      }
      
      const FinalState& ufs = apply<FinalState>(event, "UFS");
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
	  if(matched)
	    _nEtaPiPi->fill(event.weight());
	}
	else if(p.pdgId()==223) {
	  map<long,int> nRes = nCount;
	  int ncount = ntotal;
	  findChildren(p,nRes,ncount);
	  // eta pi+pi-
	  if(ncount!=1) continue;
	  bool matched = true;
	  for(auto const & val : nRes) {
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
	  if(matched)
	    _nOmegaPi->fill(event.weight());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<15;++ix) {
	unsigned int ymax=2;
       	double sigma(0.),error(0.);
       	if(ix<=4) {
       	  sigma = _nOmegaPi->val();
       	  error = _nOmegaPi->err();
       	}
	else if(ix==5) {
	  sigma = _n3Pi->val();
	  error = _n3Pi->err();
	}
	else if(ix==6) {
	  sigma = _nEtaPiPi->val();
	  error = _nEtaPiPi->err();
	}
	else if(ix==7) {
	  sigma = _n4PiC->val();
	  error = _n4PiC->err();
	}
	else if(ix==8) {
	  sigma = _n4PiN->val();
	  error = _n4PiN->err();
	}
	else if(ix==9) {
	  continue;
	}
	else if(ix==10) {
	  ymax=5;
	}
	else if(ix==11) {
	  sigma = _n2Pi->val();
	  error = _n2Pi->err();
	}
	else if(ix==12) {
	  sigma = _nKC->val();
	  error = _nKC->err();
	}
	else if(ix==13) {
	  sigma = _nKN->val();
	  error = _nKN->err();
	}
	else if(ix==14) {
	  sigma = _n5Pi->val();
	  error = _n5Pi->err();
	}
      	sigma *= crossSection()/ sumOfWeights() /nanobarn;
      	error *= crossSection()/ sumOfWeights() /nanobarn;
	for(unsigned int iy=1;iy<ymax;++iy) {
	  if(ix==10) {
	    if(iy==1) {
	      sigma = _n4PiC->val();
	      error = _n4PiC->err();
	    }
	    else if(iy==2) {
	      sigma = _n4PiN->val();
	      error = _n4PiN->err();
	    }
	    else if(iy==3) {
	      sigma = _nOmegaPi->val();
	      error = _nOmegaPi->err();
	    }
	    else {
	      sigma = _n3Pi->val();
	      error = _n3Pi->err();
	    }
	    sigma *= crossSection()/ sumOfWeights() /nanobarn;
	    error *= crossSection()/ sumOfWeights() /nanobarn;
	  }
	  Scatter2D temphisto(refData(ix, 1, iy));
	  Scatter2DPtr  mult = bookScatter2D(ix, 1, iy);
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
    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _nOmegaPi,_n2Pi,_n3Pi,_n4PiC,_n4PiN,_nEtaPiPi,_nKC,_nKN,_n5Pi;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ND_1991_I321108);


}
