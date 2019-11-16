// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMD3_2019_I1720610 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMD3_2019_I1720610);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms

      _c_all   = bookCounter("/TMP/all");
      _c_omega = bookCounter("/TMP/omega");
      _c_eta   = bookCounter("/TMP/eta");

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
      // find the final-state particles
      double weight = event.weight();
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      foreach (const Particle& p, fs.particles()) {
	nCount[p.pdgId()] += 1;
	++ntotal;
      }
      if(ntotal==7 && nCount[211]==3 && nCount[-211]==3 && nCount[111] ==1 ) {
	_c_all->fill(weight);
      }
      // find omega/phi + eta 
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      bool found=false;
      foreach (const Particle& p, ufs.particles()) {
	if(p.children().empty()) continue;
	// find the eta/omega
	if(p.pdgId()!=221 && p.pdgId()!=223) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	// eta/omega 2(pi+pi-)
	if(ncount==4) {
	  bool matched = true;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==211 ) {
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
	  if(matched) {
	    if(p.pdgId()==221)
	      _c_eta->fill(event.weight());
	    else if(p.pdgId()==223)
	      _c_omega->fill(event.weight());
	    found = true;
	    break;
	  }
	}
	if(found) break;
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ix=1;ix<4;++ix) {
	double sigma(0.),error(0.);
	if(ix==1) {
	  sigma = _c_all->val()*fact;
	  error = _c_all->err()*fact;
	}
	else if(ix==2) {
	  sigma = _c_eta->val()*fact;
	  error = _c_eta->err()*fact;
	}
	else if(ix==3) {
	  sigma = _c_omega->val()*fact;
	  error = _c_omega->err()*fact;
	}
	Scatter2D temphisto(refData(1, 1, ix));
	Scatter2DPtr  mult = bookScatter2D(1, 1, ix);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS()/MeV, x-ex2.first, x+ex2.second)) {
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
    CounterPtr _c_all,_c_omega,_c_eta;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMD3_2019_I1720610);


}
