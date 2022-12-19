// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > B*s B*s
  class CLEOIII_2006_I694869 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOIII_2006_I694869);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::abspid==533), "UFS");
      // histograms
      book(_h,1,1,1);
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for(const Particle &child : p.children()) {
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

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // extract botton hadrons
      Particles bHadrons=apply<FinalState>(event, "UFS").particles();
      for(const Particle & p1 : bHadrons) {
	if(p1.pid()<0 || p1.parents()[0].abspid()==533 ) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p1,nRes,ncount);
	bool matched=false;
	for(const Particle & p2 : bHadrons) {
	  if(p1.pid()!=-p2.pid() || p2.parents()[0].abspid()==533 ) continue;
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
	    _h->fill(10.859);
	    break;
	  }
      	}
      	if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h,crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOIII_2006_I694869);

}
