// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BELLE_2017_I1613517 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BELLE_2017_I1613517);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      _c_DpDmS   = bookCounter("/TMP/sigma_DpDmS");
      _c_DpSDmS  = bookCounter("/TMP/sigma_DpSDmS");
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
      // total hadronic and muonic cross sections
      map<long,int> nCount;
      int ntotal(0);
      foreach (const Particle& p, fs.particles()) {
	nCount[p.pdgId()] += 1;
	++ntotal;
      }
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 &&
	 ntotal==2+nCount[22])
	vetoEvent;
      // unstable charm analysis
      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      for(unsigned int ix=0;ix<ufs.particles().size();++ix) {
       	const Particle& p1 = ufs.particles()[ix];
       	int id1 = abs(p1.pdgId());
       	if(id1 != 411 && id1 != 413) continue;
      	// check fs
      	bool fs = true;
      	foreach(const Particle & child, p1.children()) {
      	  if(child.pdgId()==p1.pdgId()) {
      	    fs = false;
      	    break;
      	  }
      	}
      	if(!fs) continue;
      	// find the children
      	map<long,int> nRes = nCount;
      	int ncount = ntotal;
      	findChildren(p1,nRes,ncount);
      	bool matched=false;
       	int sign = p1.pdgId()/id1;
      	// loop over the other fs particles
      	for(unsigned int iy=ix+1;iy<ufs.particles().size();++iy) {
      	  const Particle& p2 = ufs.particles()[iy];
      	  fs = true;
      	  foreach(const Particle & child, p2.children()) {
      	    if(child.pdgId()==p2.pdgId()) {
      	      fs = false;
      	      break;
      	    }
      	  }
      	  if(!fs) continue;
       	  if(p2.pdgId()/abs(p2.pdgId())==sign) continue;
      	  int id2 = abs(p2.pdgId());
       	  if(id2 != 411 && id2 != 413) continue;
      	  if(!p2.parents().empty() && p2.parents()[0].pdgId()==p1.pdgId())
      	    continue;
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
	    if(id1==413 && id2==413) {
	      _c_DpSDmS->fill(event.weight());
	    }
	    else if((id1==411 && id2==413) ||
		    (id1==413 && id2==411)) {
	      _c_DpDmS->fill(event.weight());
	    }
	    break;
	  }
      	}
	if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights()/nanobarn;
      for(unsigned int iy=1;iy<3;++iy) {
	double sigma,error;
	if(iy==1) {
	  sigma = _c_DpDmS->val()*fact;
	  error = _c_DpDmS->err()*fact;
	}
	else if(iy==2) {
	  sigma = _c_DpSDmS->val()*fact;
	  error = _c_DpSDmS->err()*fact;
	}
	Scatter2D temphisto(refData(1, 1, iy));
        Scatter2DPtr     mult = bookScatter2D(1,1,iy);
        for (size_t b = 0; b < temphisto.numPoints(); b++) {
          const double x  = temphisto.point(b).x();
          pair<double,double> ex = temphisto.point(b).xErrs();
          pair<double,double> ex2 = ex;
          if(ex2.first ==0.) ex2. first=0.0001;
          if(ex2.second==0.) ex2.second=0.0001;
          if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
            mult   ->addPoint(x, sigma, ex, make_pair(error,error));
          }
          else {
            mult   ->addPoint(x, 0., ex, make_pair(0.,.0));
          }
        }
      }
    }
    //@}


    /// @name Histograms
    //@{
    CounterPtr _c_DpDmS, _c_DpSDmS;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BELLE_2017_I1613517);


}
