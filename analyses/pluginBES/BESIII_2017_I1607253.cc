// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+ e- > 2 phi + phi or omega
  class BESIII_2017_I1607253 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2017_I1607253);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==223 || Cuts::pid==333), "UFS");
      // counters
      for(unsigned int ix=0;ix<2;++ix)
	book(_c[ix],"TMP/c_" + toString(ix+1));
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
	if(child.children().empty()) {
	  --nRes[child.pid()];
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
      for (const Particle& p: fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      Particles phi   = ufs.particles(Cuts::pid==333);
      Particles omega = ufs.particles(Cuts::pid==223);
      if(phi.size()<2) vetoEvent;
      bool found = false;
      // first phi
      for(unsigned int ix=0;ix<phi.size();++ix) {
	const Particle & phi1 = phi[ix];
	if(phi1.children().empty()) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(phi1,nRes,ncount);
	// second phi
	for(unsigned int iy=ix+1;iy<phi.size();++iy) {
	  const Particle & phi2 = phi[iy];
	  if(phi2.children().empty()) continue;
	  map<long,int> nRes2 = nRes;
	  int ncount2 = ncount;
	  findChildren(phi2,nRes2,ncount2);
	  // search for third phi
	  for(unsigned int iz=iy+1;iz<phi.size();++iz) {
	    const Particle & phi3 = phi[iz];
	    if(phi3.children().empty()) continue;
	    map<long,int> nRes3 = nRes2;
	    int ncount3 = ncount2;
	    findChildren(phi3,nRes3,ncount3);
	    found = true;
	    for(auto const & val : nRes3) {
	      if(val.second!=0) {
		found = false;
		break;
	      }
	    }
	    if(found) {
	      _c[1]->fill();
	      break;
	    }
	  }
	  if(found) break;
	  // search for omega
	  for(unsigned int iz=0;iz<omega.size();++iz) {
	    const Particle & omega1 = omega[iz];
	    if(omega1.children().empty()) continue;
	    map<long,int> nRes3 = nRes2;
	    int ncount3 = ncount2;
	    findChildren(omega1,nRes3,ncount3);
	    found = true;
	    for(auto const & val : nRes3) {
	      if(val.second!=0) {
		found = false;
		break;
	      }
	    }
	    if(found) {
	      _c[0]->fill();
	      break;
	    }
	  }
	  if(found) break;
	}
	if(found) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact =  crossSection()/ sumOfWeights() /femtobarn;
      for(unsigned int ix=0;ix<2;++ix) {
	double sigma = _c[ix]->val()*fact;
	double error = _c[ix]->err()*fact;
	Scatter2D temphisto(refData(1+ix, 1, 1));
	Scatter2DPtr  mult;
	book(mult, 1+ix, 1, 1);
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

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2017_I1607253);

}
