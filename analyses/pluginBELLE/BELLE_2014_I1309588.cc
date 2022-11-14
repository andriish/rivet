// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > pi+ pi- pi0 chi_b
  class BELLE_2014_I1309588 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2014_I1309588);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<3;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h[ix][iy],1+ix,1,1+iy);
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
      // loop over any chi mesons
      for(const Particle & chi : ufs.particles(Cuts::pid==20553 or
					       Cuts::pid==555)) {
	if(chi.children().empty()) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(chi,nRes,ncount);
	// first check for 3 pi chi_b
	bool matched3Pi = false;
	if(ncount==3) {
	  for(auto const & val : nRes) {
	    if(abs(val.first)==PID::PIPLUS || val.first==PID::PI0) {
	      if (val.second!=1) {
		matched3Pi = false;
		break;
	      }
	    }
	    else if (val.second!=0) {
	      matched3Pi = false;
	      break;
	    }
	  }
	}
	// then for omega chi_b
	bool matchedOmega = false;
      	for(const Particle & omega : ufs.particles(Cuts::pid==223)) {
       	  map<long,int> nRes2 = nRes;
       	  int ncount2 = ncount;
       	  findChildren(omega,nRes2,ncount2);
       	  matchedOmega = true;
       	  for(auto const & val : nRes2) {
       	    if (val.second!=0) {
       	      matchedOmega = false;
       	      break;
       	    }
       	  }
	  if (matchedOmega) break;
	}
	if(!matched3Pi && !matchedOmega) continue;
	unsigned int iloc= chi.pid()==20443 ? 0 : 1;
	if(matched3Pi)
	  _h[0][iloc]->fill(10.867);
	if(matchedOmega)
	  _h[1][iloc]->fill(10.867);
	if(matched3Pi && !matchedOmega)
	  _h[2][iloc]->fill(10.867);
	break;
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      double fact =  crossSection()/ sumOfWeights() /picobarn;
      for(unsigned int ix=0;ix<3;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  scale(_h[ix][iy],fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2014_I1309588);

}
