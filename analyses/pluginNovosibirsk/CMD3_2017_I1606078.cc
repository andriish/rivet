// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief CMD3 pi+pi-pi0eta cross section
  class CMD3_2017_I1606078 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMD3_2017_I1606078);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");


      _c_all   = bookCounter("/TMP/all");
      _c_omega = bookCounter("/TMP/omega");
      _c_rho   = bookCounter("/TMP/rho");
      _c_other   = bookCounter("/TMP/other");
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
      // find omega/phi + eta 
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      bool found = false, foundOmegaPhi = false;
      foreach (const Particle& p, ufs.particles()) {
	if(p.children().empty()) continue;
	// find the eta
	if(p.pdgId()!=221) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	// eta pi+pi-pi0
	if(ncount==3) {
	  bool matched = true;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==211 || val.first==111 ) {
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
	    _c_all->fill(event.weight());
	    found = true;
	  }
	}
	foreach (const Particle& p2, ufs.particles()) {
	  if(p2.pdgId()!=223 && p2.pdgId()!=333) continue;
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
	    if(p2.pdgId()==223)
	      _c_omega->fill(event.weight());
	    foundOmegaPhi=true;
	  }
	}
      }
      // find a_0 rho
      bool founda0Rho=false;
      foreach (const Particle& p, ufs.particles()) {
      	if(p.children().empty()) continue;
	// find the rho
	if(p.pdgId()!=113 && abs(p.pdgId())!=213)
	  continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
	int a1id(0);
	if(p.pdgId()==213)
	  a1id = 9000211;
	else if(p.pdgId()==-213)
	  a1id =-9000211;
	else
	  a1id = 9000111;
	foreach (const Particle& p2, ufs.particles()) {
	  if(p2.pdgId()!=a1id) continue;
	  map<long,int> nResB = nRes;
	  int ncountB = ncount;
	  findChildren(p2,nResB,ncountB);
	  if(ncountB!=0) continue;
	  bool matched = true;
	  for(auto const & val : nResB) {
	    if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) {
	    _c_rho->fill(event.weight());
	    founda0Rho=true;
	    break;
	  }
	}
	if(founda0Rho) break;
      }
      if(found && !foundOmegaPhi && ! founda0Rho)
	_c_other->fill(weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ix=1;ix<5;++ix) {
	double sigma(0.),error(0.);
	if(ix==1) {
	  sigma = _c_all->val()*fact;
	  error = _c_all->err()*fact;
	}
	else if(ix==2) {
	  sigma = _c_omega->val()*fact;
	  error = _c_omega->err()*fact;
	}
	else if(ix==3) {
	  sigma = _c_rho->val()*fact;
	  error = _c_rho->err()*fact;
	}
	else if(ix==4) {
	  sigma = _c_other->val()*fact;
	  error = _c_other->err()*fact;
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
    CounterPtr _c_all, _c_omega, _c_rho, _c_other;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMD3_2017_I1606078);


}
