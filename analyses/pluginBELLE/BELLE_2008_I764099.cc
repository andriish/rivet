// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BELLE_2008_I764099 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BELLE_2008_I764099);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      _nUps1pipi = bookCounter("TMP/nUps1pipi");
      _nUps2pipi = bookCounter("TMP/nUps2pipi");
      _nUps3pipi = bookCounter("TMP/nUps3pipi");
      _nUps1KK   = bookCounter("TMP/nUps1KK");
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
      foreach (const Particle& p, ufs.particles()) {
	if(p.children().empty()) continue;
	if(p.pdgId() !=   553  &&
	   p.pdgId() != 100553 &&
	   p.pdgId() != 200553 ) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	if(ncount!=2) continue;
	bool matched = true;
	for(auto const & val : nRes) {
	  if(abs(val.first)==321 || abs(val.first)==211) {
	    continue;
	  }
	  else if(val.second!=0) {
	    matched = false;
	    break;
	  }
	}
	if(matched) {
	  if(nRes[211]==1 && nRes[-211]==1 ) {
	    if(p.pdgId()==553)
	      _nUps1pipi->fill(event.weight());
	    if(p.pdgId()==100553)
	      _nUps2pipi->fill(event.weight());
	    if(p.pdgId()==200553)
	      _nUps3pipi->fill(event.weight());
	  }
	  else if(nRes[321]==1 && nRes[-321]==1) {
	    if(p.pdgId()==553)
	      _nUps1KK->fill(event.weight());
	  }	  
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /picobarn;
      for(unsigned int ix=1;ix<5;++ix) {
	double sigma,error;
	if(ix==1) {
	  sigma = _nUps1pipi->val()*fact;
	  error = _nUps1pipi->err()*fact;
	}
	else if(ix==2) {
	  sigma = _nUps2pipi->val()*fact;
	  error = _nUps2pipi->err()*fact;
	}
	else if(ix==3) {
	  sigma = _nUps3pipi->val()*fact;
	  error = _nUps3pipi->err()*fact;
	}
	else if(ix==4) {
	  sigma = _nUps1KK->val()*fact;
	  error = _nUps1KK->err()*fact;
	}
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
    CounterPtr _nUps1pipi,_nUps2pipi,_nUps3pipi,_nUps1KK;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BELLE_2008_I764099);

}
