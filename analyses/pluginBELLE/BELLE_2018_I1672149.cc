// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief gamma gamma -> eta' pi+pi-
  class BELLE_2018_I1672149 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2018_I1672149);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // counters
      book(_nCount[0],"TMP/nf2EtaPrime");
      book(_nCount[1],"TMP/nEtaPrimePiPi");
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
      // first check for f_2 eta'
      bool foundf2(false);
      Particles f2s = ufs.particles(Cuts::pid==225);
      Particles etaps = ufs.particles(Cuts::pid==331);
      for(const Particle & f2 : f2s) {
	bool matched=false;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(f2,nRes,ncount);
	for(const Particle & etap : etaps) {
	  map<long,int> nRes2=nRes;
	  int ncount2 = ncount;
	  findChildren(etap,nRes2,ncount2);
	  if(ncount2 !=0 ) continue;
	  matched = true;
	  for(auto const & val : nRes2) {
	    if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) {
	    foundf2 = true;
	    break;
	  }
	}
	if (foundf2) break;
      }
      // if we have the f_2 eta' state
      if(foundf2) {
	_nCount[0]->fill();
	if(sqrtS()>2.26*GeV) vetoEvent;
      }
      // see if we have eta' pi+pi-
      bool foundetap(false);
      for(const Particle & etap : etaps) {
	bool matched=false;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(etap,nRes,ncount);
	if(ncount !=2 ) continue;
	matched=true;
	for(auto const & val : nRes) {
	  if(abs(val.first)==211) {
	    if(val.second!=1) {
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
	  foundetap = true;
	  break;
	}
      }
      // check there's no eta_c
      if(foundetap && sqrtS()>2.62 && sqrtS()<3.06) {
	for(const Particle & etac : ufs.particles(Cuts::pid==441)) {
	  bool matched=false;
	  map<long,int> nRes=nCount;
	  int ncount = ntotal;
	  findChildren(etac,nRes,ncount);
	  for(auto const & val : nRes) {
	    if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) {
	    foundetap=false;
	    break;
	  }
	}
      }
      if(foundetap) _nCount[1]->fill();
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      // loop over tables in paper
      for(unsigned int ix=0;ix<2;++ix) {
	double sigma = _nCount[ix]->val()*fact;
	double error = _nCount[ix]->err()*fact;
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
    CounterPtr _nCount[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2018_I1672149);

}
