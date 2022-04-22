// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief gamma gamma -> K*0K*0
  class ARGUS_1987_I248680 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1987_I248680);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // book histos
      if(inRange(sqrtS()/GeV,1.6,3.5)) {
	for(unsigned int ix=0;ix<9;++ix)
	  book(_nMeson[ix],"TMP/nMeson_"+toString(ix+1));
      }
      else
	throw Error("Invalid CMS energy for ARGUS_1987_I248680");
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
      // find any K* mesons
      int ires=-1;
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles Kstar=ufs.particles(Cuts::abspid==313);
      for (unsigned int ix=0;ix<Kstar.size();++ix) {
       	if(Kstar[ix].children().empty()) continue;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(Kstar[ix],nRes,ncount);
	bool matched=false;
	// K*K*
	for (unsigned int iy=ix+1;iy<Kstar.size();++iy) {
	  if(Kstar[iy].children().empty()) continue;
	  if(Kstar[ix].pid()!=-Kstar[iy].pid()) continue;
	  map<long,int> nRes2=nRes;
	  int ncount2 = ncount;
	  findChildren(Kstar[iy],nRes2,ncount2);
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
	  _nMeson[1]->fill();
	  ires=7;
	  break;
	}
	int sign = Kstar[ix].pid()/Kstar[ix].abspid();
	// three body intermediate states
	if(ncount==2) {
	  // K*0 K- pi+ +ccd
	  matched=true;
	  for(auto const & val : nRes) {
	    if (val.first==sign*211 || val.first==-sign*321) {
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
	    _nMeson[2]->fill();
	    ires=6;
	    break;
	  }
	}
      }
      // look for phi modes
      for (const Particle & p : ufs.particles(Cuts::pid==PID::PHI)) {
       	if(p.children().empty()) continue;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(p,nRes,ncount);
	if(ncount==2) {
	  bool matched=true;
	  for(auto const & val : nRes) {
	    if (abs(val.first)==211) {
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
	    ires=8;
	    break;
	  }
	}
      }
      // 4 meson modes
      if(ntotal==4 &&
	 nCount[PID::KPLUS ]==1 && nCount[PID::KMINUS ]==1 &&
	 nCount[PID::PIPLUS]==1 && nCount[PID::PIMINUS]==1 ) {
	_nMeson[0]->fill();
	_nMeson[4]->fill();
	if(ires<0) {
	  _nMeson[3]->fill();
	  _nMeson[5]->fill();
	}
	else _nMeson[ires]->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      // loop over tables in paper
      for(unsigned int ih=1;ih<6;++ih) {
       	unsigned int imax= ih!=5 ? 2 : 6;
       	for(unsigned int iy=1;iy<imax;++iy) {
	  unsigned int iloc=ih+iy-2;
       	  assert(iloc<=8);
       	  double sigma = _nMeson[iloc]->val()*fact;
       	  double error = _nMeson[iloc]->err()*fact;
       	  Scatter2D temphisto(refData(ih, 1, iy));
       	  Scatter2DPtr mult;
       	  book(mult, ih, 1, iy);
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
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nMeson[9];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1987_I248680);

}
