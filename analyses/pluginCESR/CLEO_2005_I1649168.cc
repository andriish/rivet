// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CLEO_2005_I1649168 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2005_I1649168);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_4S, 1, 1, 7);
      book(_h_5S, 1, 1, 8);
      book(_c_4S, "/TMP/N4S");
      book(_c_5S, "/TMP/N5S");
    }


    void findDecayProducts(const Particle & p, Particles & ds) {
      for(const Particle & child : p.children()) {
	if(child.abspid()==431)
	  ds.push_back(child);
	else if(!child.children().empty())
	  findDecayProducts(child,ds);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the Upsilons among the unstables
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==300553 || Cuts::pid==400553);
      for (const Particle& ups : upsilons) {
        LorentzTransform cms_boost;
        if (ups.p3().mod() > 1*MeV)
          cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
	Particles ds;
        findDecayProducts(ups, ds);

	if(ups.pid()==300553)
	  _c_4S->fill();
	else
	  _c_5S->fill();

	for(const Particle & p : ds) {
          FourMomentum p2 = cms_boost.transform(p.momentum());
	  double x = 2.*p2.p3().mod()/ups.mass();
	if(ups.pid()==300553)
	  _h_4S->fill(x);
	else
	  _h_5S->fill(x);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      if(_h_4S->effNumEntries()!=0)
	scale(_h_4S   ,100./ *_c_4S);
      if(_h_5S->effNumEntries()!=0)
	scale(_h_5S   ,100./ *_c_5S);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_4S,_h_5S;
    CounterPtr _c_4S,_c_5S;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CLEO_2005_I1649168);


}
