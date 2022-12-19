// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"

namespace Rivet {


  /// @brief e+ e- > D Dbar
  class CLEOC_2007_I761696 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOC_2007_I761696);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::abspid==411 ||
				Cuts::abspid==421), "UFS");
      // histograms
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1,1,1+ix);
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

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for(const Particle & p1 : ufs.particles()) {
	if(p1.pid()<0) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p1,nRes,ncount);
	bool matched=false;
	for(const Particle & p2 : ufs.particles()) {
     	  if(p2.pid()!=-p1.pid()) continue;
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
   	  if(matched) break;
   	}
   	if(matched) {
	  _h[2]->fill(sqrtS()/MeV);
	  if(p1.abspid()==421)      _h[0]->fill(sqrtS()/MeV);
	  else if(p1.abspid()==411) _h[1]->fill(sqrtS()/MeV);
	  break;
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/sumOfWeights()/nanobarn;
      Scatter2DPtr ratio;
      book(ratio,1,1,4);
      divide(_h[1],_h[0],ratio);
      for(unsigned int ix=0;ix<3;++ix) scale(_h[ix],fact*2.);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOC_2007_I761696);

}
