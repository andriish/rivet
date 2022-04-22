// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief gamma gamma -> omega omega
  class ARGUS_1996_I403304 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1996_I403304);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // book histos
      if(inRange(sqrtS()/GeV,1.6,3.0)) {
	for(unsigned int ix=0;ix<3;++ix)
	  book(_nMeson[ix],"TMP/nMeson_"+toString(ix+1));
      }
      else
	throw Error("Invalid CMS energy for ARGUS_1996_I403304");
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
      // find omega mesons
      Particles omega=apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==223);
      bool nonRes=true;
      for (unsigned int ix=0;ix<omega.size();++ix) {
       	if(omega[ix].children().empty()) continue;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(omega[ix],nRes,ncount);
	bool matched=false;
	// omega omega
	for (unsigned int iy=ix+1;iy<omega.size();++iy) {
	  if(omega[iy].children().empty()) continue;
	  map<long,int> nRes2=nRes;
	  int ncount2 = ncount;
	  findChildren(omega[iy],nRes2,ncount2);
	  if(ncount2 !=0 ) continue;
	  matched = true;
	  for(auto const & val : nRes2) {
	    if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched)
	    break;
	}
	if(matched) {
	  _nMeson[0]->fill();
	  nonRes=false;
	  break;
	}
	matched=true;
	for(auto const & val : nRes) {
	  if (abs(val.first)==211 || val.first==111) {
	    if(val.second!=1) {
	      matched = false;
	      break;
	    }
	  }
	  else {
	    if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	}
	if(matched) {
	  nonRes=false;
	  _nMeson[1]->fill();
	  break;
	}
      }
      if(nonRes && ntotal==6 && nCount[PID::PI0]==2 &&
	 nCount[PID::PIPLUS]==2 && nCount[PID::PIMINUS]==2 ) {
	_nMeson[2]->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      // loop over tables in paper
      for(unsigned int ix=0;ix<3;++ix) {
	double sigma = _nMeson[ix]->val()*fact;
	double error = _nMeson[ix]->err()*fact;
	Scatter2D temphisto(refData(1, 1, 1+ix));
	Scatter2DPtr mult;
	book(mult, 1, 1, 1+ix);
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
    CounterPtr _nMeson[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1996_I403304);

}
