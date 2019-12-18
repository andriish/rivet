// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D*2 spectrum and decay angle
  class ARGUS_1989_I268577 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ARGUS_1989_I268577);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_D2_x     , 2, 1, 1);
      book(_h_D2_ctheta, 3, 1, 1);
    }

   
   /// Recursively walk the decay tree to find decay products of @a p
   void findDecayProducts(Particle mother, Particles & dstar, Particles & pi,unsigned int & ncount) {
     for(const Particle & p: mother.children()) {
	if(p.abspid()==413)
	  dstar.push_back(p);
	else if(p.abspid()==211)
	  pi.push_back(p);
	ncount +=1;
     }
   }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==425)) {
	const double xp = 2.*p.p3().mod()/sqrtS();
	_h_D2_x->fill(xp);
	// decay products
	Particles dstar,pi;
	unsigned int ncount=0;
	findDecayProducts(p,dstar,pi,ncount);
	if(ncount!=2 || dstar.size()!=1 || pi.size()!=1) continue;
	if(dstar[0].pid()/p.pid()<0) continue;
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	Vector3 axis = boost.transform(dstar[0].momentum()).p3().unit();
	double cosL  = axis.dot(p.momentum().p3().unit());
	// decay angles
	_h_D2_ctheta->fill(cosL);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_D2_x);
      normalize(_h_D2_ctheta);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_D2_x,_h_D2_ctheta;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ARGUS_1989_I268577);


}
