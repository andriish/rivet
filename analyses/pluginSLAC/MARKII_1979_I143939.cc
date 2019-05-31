// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MARKII_1979_I143939 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MARKII_1979_I143939);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      _c_hadrons = bookCounter("/TMP/sigma_hadrons");
      _c_muons   = bookCounter("/TMP/sigma_muons");
      _c_DD      = bookCounter("/TMP/sigma_DD");
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
	_c_muons->fill(event.weight());
      // everything else
      else
	_c_hadrons->fill(event.weight());
      // identified final state with D mesons
      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      for(unsigned int ix=0;ix<ufs.particles().size();++ix) {
	bool matched = false;
       	const Particle& p1 = ufs.particles()[ix];
       	int id1 = abs(p1.pdgId());
       	if(id1 != 411 && id1 != 421) continue;
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
       	  if(p2.pdgId()/abs(p2.pdgId())==p1.pdgId()/abs(p1.pdgId())) continue;
       	  int id2 = abs(p2.pdgId());
       	  if(id2 != 411 && id2 != 421) continue;
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
	    _c_DD  ->fill(event.weight());
	    break;
	  }
	}
	if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      Scatter1D R = *_c_hadrons/ *_c_muons;
      double              rval = R.point(0).x();
      pair<double,double> rerr = R.point(0).xErrs();
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      double sig_h = _c_hadrons->val()*fact;
      double err_h = _c_hadrons->err()*fact;
      double sig_m = _c_muons  ->val()*fact;
      double err_m = _c_muons  ->err()*fact;
      Scatter2D temphisto(refData(2, 1, 1));
      Scatter2DPtr hadrons  = bookScatter2D("sigma_hadrons");
      Scatter2DPtr muons    = bookScatter2D("sigma_muons"  );
      Scatter2DPtr     mult = bookScatter2D(2,1,1);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
	const double x  = temphisto.point(b).x();
	pair<double,double> ex = temphisto.point(b).xErrs();
	pair<double,double> ex2 = ex;
	if(ex2.first ==0.) ex2. first=0.0001;
	if(ex2.second==0.) ex2.second=0.0001;
	if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
	  mult   ->addPoint(x, rval, ex, rerr);
	  hadrons->addPoint(x, sig_h, ex, make_pair(err_h,err_h));
	  muons  ->addPoint(x, sig_m, ex, make_pair(err_m,err_m));
	}
	else {
	  mult   ->addPoint(x, 0., ex, make_pair(0.,.0));
	  hadrons->addPoint(x, 0., ex, make_pair(0.,.0));
	  muons  ->addPoint(x, 0., ex, make_pair(0.,.0));
	}
      }
      double sigma = _c_DD->val()*fact;
      double error = _c_DD->err()*fact;
      Scatter2D temphisto2(refData(3, 1, 1));
      mult = bookScatter2D(3,1,1);
      for (size_t b = 0; b < temphisto2.numPoints(); b++) {
	const double x  = temphisto2.point(b).x();
	pair<double,double> ex = temphisto2.point(b).xErrs();
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

    //@}


    /// @name Histograms
    //@{
    CounterPtr _c_hadrons, _c_muons, _c_DD;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MARKII_1979_I143939);


}
