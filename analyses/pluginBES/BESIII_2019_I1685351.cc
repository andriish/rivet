// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > mu+ mu-, pi+ pi- eta and 5pi near J/psi
  class BESIII_2019_I1685351 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1685351);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==PID::ETA), "UFS");
      // histograms
      for(unsigned int ix=0;ix<3;++ix)
	book(_c[ix],"TMP/c_"+toString(ix+1));
    }
    
    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for(const Particle &child : p.children()) {
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
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // mu+mu- + photons
      if(nCount[-13]==1 && nCount[13]==1 && ntotal==2+nCount[22])
	_c[0]->fill();
      else if (nCount[111]==1 && nCount[211]==2 && nCount[-211]==2  && ntotal==5+nCount[22])
	_c[2]->fill();
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      // loop over eta mesons
      for (const Particle& p : ufs.particles()) {
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	bool matched = true;
	for(auto const & val : nRes) {
	  if(abs(val.first)==211) {
	    if(val.second !=1) {
	      matched = false;
	      break;
	    }
	  }
	  else if(val.first!=PID::PHOTON && val.second!=0) {
	    matched = false;
	    break;
	  }
	}
	if(!matched) continue;
	_c[1]->fill();
	break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      for(unsigned int ix=0;ix<3;++ix) {
	double sig_m = _c[ix]  ->val()*fact;
	double err_m = _c[ix]  ->err()*fact;
	Scatter2D temphisto(refData(1, 1, 1+ix));
	Scatter2DPtr cross;
	book(cross, 1, 1, 1+ix);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS()/MeV, x-ex2.first, x+ex2.second))
	    cross  ->addPoint(x, sig_m, ex, make_pair(err_m,err_m));
	  else
	    cross  ->addPoint(x, 0., ex, make_pair(0.,.0));
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1685351);

}