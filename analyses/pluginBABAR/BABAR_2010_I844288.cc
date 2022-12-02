// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief gamma gamma -> D Dbar
  class BABAR_2010_I844288 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2010_I844288);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<double>("PID", 100445);
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(),"FS");
      declare(UnstableParticles(Cuts::abspid==411 ||
				Cuts::abspid==421), "UFS");
      // histograms
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1+ix,1,1);
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if (child.children().empty()) {
          --nRes[child.pid()];
          --ncount;
        } else {
          findChildren(child,nRes,ncount);
        }
      }
    }

    bool findScattered(Particle beam, double& q2) {
      bool found = false;
      Particle scat = beam;
      while (!scat.children().empty()) {
        found = false;
        for (const Particle & p : scat.children()) {
          if (p.pid()==scat.pid()) {
            scat=p;
            found=true;
            break;
          }
        }
        if (!found) break;
      }
      if (!found) return false;
      q2 = -(beam.momentum() - scat.momentum()).mass2();
      return true;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find scattered leptons and calc Q2
      const Beam& beams = apply<Beam>(event, "Beams");
      double q12 = -1, q22 = -1;
      if (!findScattered(beams.beams().first,  q12)) vetoEvent;
      if (!findScattered(beams.beams().second, q22)) vetoEvent;
      // check the final state
      const FinalState & fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      // find the meson
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p1 : ufs.particles()) {
	if(p1.children().empty() || p1.pid()<0) continue;
	bool matched=false;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p1,nRes,ncount);
	for (const Particle& p2 : ufs.particles(Cuts::pid==-p1.pid())) {
	  map<long,int> nRes2 = nRes;
	  int ncount2 = ncount;
	  findChildren(p2,nRes2,ncount2);
	  matched = true;
	  for(auto const & val : nRes2) {
	    if(abs(val.first)==11) {
	      if(val.second!=1) {
		matched = false;
		break;
	      }
	    }
	    else if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if (matched) {
	    FourMomentum pDD = p1.momentum()+p2.momentum();
	    _h[0]->fill(pDD.mass());
	    if(p1.parents()[0].pid()==_pid && p2.parents()[0].pid()==_pid) {
	      Vector3 axis = pDD.p3().unit();
	      LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(pDD.betaVec());
	      double cTheta  = abs(axis.dot(boost.transform(p1).p3().unit()));
	      _h[1]->fill(cTheta);
	    }
	    break;
	  }
	}
	if(matched) break;
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    int _pid;
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2010_I844288);

}
