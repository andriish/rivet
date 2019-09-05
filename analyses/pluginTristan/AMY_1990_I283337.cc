// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Hemispheres.hh"

namespace Rivet {

  /// @brief Event shapes at 55.2
  class AMY_1990_I283337 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(AMY_1990_I283337);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      const FinalState fs;
      declare(fs, "FS");
      const Sphericity sphere(fs);
      declare(sphere, "Sphericity");
      const Thrust thrust(fs);
      declare(thrust, "Thrust");
      declare(Hemispheres(sphere), "Hemispheres");
      // histograms
      _histRapidityT  = bookHisto1D( 1, 1, 1);
      _histScaledMom  = bookHisto1D( 2, 1, 1);
      _histPl         = bookHisto1D( 3, 1, 1);
      _histPt         = bookHisto1D( 4, 1, 1);
      _histPt2        = bookHisto1D( 5, 1, 1);
      _histPtIn       = bookHisto1D( 6, 1, 1);
      _histPtOut      = bookHisto1D( 7, 1, 1);
      _histMeanPtIn2  = bookHisto1D( 8, 1, 1);
      _histMeanPtOut2 = bookHisto1D( 9, 1, 1);
      _histNtheta     = bookHisto1D(10, 1, 1);
      _histEtheta     = bookHisto1D(11, 1, 1);
      _histThrust     = bookHisto1D(12, 1, 1);
      _histMajor      = bookHisto1D(13, 1, 1);
      _histMinor      = bookHisto1D(14, 1, 1);
      _histOblateness = bookHisto1D(15, 1, 1);
      _histSphericity = bookHisto1D(16, 1, 1);
      _histAplanarity = bookHisto1D(17, 1, 1);
      _histQx         = bookHisto1D(18, 1, 1);
      _histQ21        = bookHisto1D(19, 1, 1);
      _histRhoLight   = bookHisto1D(20, 1, 1);
      _histRhoHeavy   = bookHisto1D(21, 1, 1);
      _histRhoDiff    = bookHisto1D(22, 1, 1);
      _wSum=0.;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");
      const double weight = event.weight();
      _wSum += weight;

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      // Thrusts
      MSG_DEBUG("Calculating thrust");
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      _histThrust    ->fill(thrust.thrust()     , weight);
      _histMajor     ->fill(thrust.thrustMajor(), weight);
      _histMinor     ->fill(thrust.thrustMinor(), weight);
      _histOblateness->fill(thrust.oblateness() , weight);
      // Sphericities
      MSG_DEBUG("Calculating sphericity");
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      _histSphericity->fill(sphericity.sphericity(), weight);
      _histAplanarity->fill(sphericity.aplanarity(), weight);
      _histQx        ->fill(sqrt(1./3.)*(sphericity.lambda1()-sphericity.lambda2()), weight);
      _histQ21       ->fill(sphericity.lambda2()-sphericity.lambda3(), weight);
      // Hemispheres
      MSG_DEBUG("Calculating hemisphere variables");
      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");
      _histRhoHeavy->fill(hemi.scaledM2high(), weight);
      _histRhoLight->fill(hemi.scaledM2low() , weight);
      _histRhoDiff ->fill(hemi.scaledM2diff(), weight);
      // single particle distributions
      double pTIn2(0.),pTOut2(0.);
      unsigned int nCharged(0);
      for (const Particle& p : fs.particles()) {
        // Get momentum and energy of each particle.
        const Vector3 mom3 = p.p3();
        const double energy = p.E();
        const double mom = mom3.mod();
        const double scaledMom = mom/meanBeamMom;
        const double momT = dot(thrust.thrustAxis(), mom3);
        const double momS = dot(sphericity.sphericityAxis(), mom3);
        const double pTinS = dot(mom3, sphericity.sphericityMajorAxis());
        const double pToutS = dot(mom3, sphericity.sphericityMinorAxis());
        const double pT = sqrt(pow(pTinS, 2) + pow(pToutS, 2));

        const double rapidityT = 0.5 * std::log((energy + momT) / (energy - momT));
	double angle = sphericity.sphericityAxis().angle(p.p3())/M_PI*180.;
	if(angle>90.) angle=180.-angle;
	if(PID::isCharged(p.pdgId())) {
	  _histScaledMom->fill(scaledMom, weight);
	  _histRapidityT->fill(fabs(rapidityT), weight);
	  _histPl       ->fill(fabs(momS)     , weight);
	  _histPt       ->fill(pT             , weight);
	  _histPt2      ->fill(sqr(pT)        , weight);
	  _histPtIn     ->fill(fabs(pTinS)    , weight);
	  _histPtOut    ->fill(fabs(pToutS)   , weight);
	  pTIn2  += sqr(pTinS);
	  pTOut2 += sqr(pToutS);
	  _histNtheta->fill(angle,weight);
	  ++nCharged;
	}
	_histEtheta->fill(angle,energy*weight); 
      }
      _histMeanPtIn2 ->fill( pTIn2/nCharged,weight);
      _histMeanPtOut2->fill(pTOut2/nCharged,weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // histograms
      scale(_histRapidityT , 1./_wSum);
      scale(_histScaledMom , 1./_wSum);
      scale(_histPl        , 1./_wSum);
      scale(_histPt        , 1./_wSum);
      scale(_histPt2       , 1./_wSum);
      scale(_histPtIn      , 1./_wSum);
      scale(_histPtOut     , 1./_wSum);
      scale(_histMeanPtIn2 , 1./_wSum);
      scale(_histMeanPtOut2, 1./_wSum);
      scale(_histNtheta    , 1./_wSum);
      scale(_histEtheta    , 1./_wSum);
      scale(_histThrust    , 1./_wSum);
      scale(_histMajor     , 1./_wSum);
      scale(_histMinor     , 1./_wSum);
      scale(_histOblateness, 1./_wSum);
      scale(_histSphericity, 1./_wSum);
      scale(_histAplanarity, 1./_wSum);
      scale(_histQx        , 1./_wSum);
      scale(_histQ21       , 1./_wSum);
      scale(_histRhoLight  , 1./_wSum);
      scale(_histRhoHeavy  , 1./_wSum);
      scale(_histRhoDiff   , 1./_wSum);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _histRapidityT, _histScaledMom, _histPl, _histPt, _histPt2, _histPtIn, _histPtOut,
      _histMeanPtIn2, _histMeanPtOut2, _histNtheta, _histEtheta, _histThrust, _histMajor, _histMinor,
      _histOblateness, _histSphericity, _histAplanarity, _histQx, _histQ21, _histRhoLight,
      _histRhoHeavy, _histRhoDiff;
    double _wSum;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(AMY_1990_I283337);


}
