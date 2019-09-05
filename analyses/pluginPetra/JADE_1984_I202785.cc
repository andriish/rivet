// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief D*+ production at 34.4 GeV
  class JADE_1984_I202785 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(JADE_1984_I202785);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      _h_x     = bookHisto1D(1, 1, 1);
      _h_theta = bookHisto1D(3, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get event weight for histo filling
      const double weight = event.weight();
      
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      Vector3 axis = beams.first.pdgId() == PID::EMINUS ?
	beams.first.p3().unit() : beams.second.p3().unit();
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      foreach (const Particle& p, ufs.particles(Cuts::abspid==413)) {
	double modp = p.p3().mod();
	double xP =modp/meanBeamMom;
	_h_x->fill(xP, weight);
	if(xP>0.4&&p.pdgId()>0) {
	  _h_theta->fill(p.p3().dot(axis)/modp,weight);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_theta,33.6,false); // normalize to data
      scale(_h_x, crossSection()/microbarn/sumOfWeights()*sqr(sqrtS()));

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_x, _h_theta;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(JADE_1984_I202785);


}
