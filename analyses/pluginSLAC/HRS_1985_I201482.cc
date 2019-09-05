// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {


  /// @brief event shapes at 29 GeV
  class HRS_1985_I201482 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(HRS_1985_I201482);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      const ChargedFinalState cfs;
      declare(cfs, "FS");
      declare(Sphericity(cfs), "Sphericity");
      const Thrust thrust(cfs);
      declare(thrust, "Thrust");

      // Book histograms
      _histSphericity = bookHisto1D( 1, 1, 1);
      _histThrust     = bookHisto1D( 3, 1, 1);
      _histThrust2Jet = bookHisto1D( 4, 1, 1);
      _histAplanarity = bookHisto1D( 6, 1, 1);
      _histZ          = bookHisto1D(10, 1, 1);
      _histZ2Jet      = bookHisto1D(11, 1, 1);
      _histZScale     = bookHisto1D(12, 1, 1);
      _histZJet[0]    = bookHisto1D(13, 1, 1);
      _histZJet[1]    = bookHisto1D(14, 1, 1);
      _histZJet[2]    = bookHisto1D(15, 1, 1);
      _histXFeyn      = bookHisto1D(16, 1, 1);
      _histXFeyn2Jet  = bookHisto1D(17, 1, 1);
      _histRap        = bookHisto1D(19, 1, 1);
      _histRap2Jet    = bookHisto1D(20, 1, 1);
      _histPtT        = bookHisto1D(22, 1, 1);
      _histPtT2Jet    = bookHisto1D(23, 1, 1);
      _histPtTIn      = bookHisto1D(24, 1, 1);
      _histPtTOut     = bookHisto1D(25, 1, 1);
      _wSum=0.;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // event weight
      const double weight = event.weight();
      // require 5 charged particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();
      if(numParticles<5) vetoEvent;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      // calc thrust and sphericity
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      Vector3 axis = thrust.thrustAxis();
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      // identify two and three jet regions
      bool twoJet   = sphericity.sphericity()<=0.25 && sphericity.aplanarity()<=0.1;
      //bool threeJet = sphericity.sphericity() >0.25 && sphericity.aplanarity()<=0.1;
      _wSum  += weight;
      if(twoJet) _wSum2 += weight;
      // basic event shapes
      _histSphericity->fill(sphericity.sphericity(),weight);
      _histThrust    ->fill(thrust.thrust(),weight);
      _histAplanarity->fill(sphericity.aplanarity(),weight);
      if(twoJet)
	_histThrust2Jet->fill(thrust.thrust(),weight);
      double pTSqIn  = 0.;
      double pTSqOut = 0.;
      unsigned int iPlus(0),iMinus(0);
      // single particle  dists
      for(const Particle & p : sortBy(fs.particles(),cmpMomByP)) {
	const double z  = p.p3().mod()/meanBeamMom;
	const double momT = axis.dot(p.p3());
	const double xF = fabs(momT)/meanBeamMom;
        const double energy = p.E();
        const double rap = 0.5 * std::log((energy + momT) / (energy - momT));
        const double pTin  = dot(p.p3(), thrust.thrustMajorAxis());
        const double pTout = dot(p.p3(), thrust.thrustMinorAxis());
	const double pT2 = sqr(pTin)+sqr(pTout);
	pTSqIn  += sqr(dot(p.p3(), sphericity.sphericityMajorAxis()));
	pTSqOut += sqr(dot(p.p3(), sphericity.sphericityMinorAxis()));
	_histZ     ->fill(z         ,  weight);
	_histZScale->fill(z         ,  weight);
	_histXFeyn ->fill(xF        ,z*weight);
	_histRap   ->fill(rap       ,  weight);
	_histPtT   ->fill(pT2       ,  weight);
	if(twoJet) {
	  _histZ2Jet    ->fill(z  ,  weight);
	  _histXFeyn2Jet->fill(xF ,z*weight);
	  _histRap2Jet  ->fill(rap,  weight);
	  _histPtT2Jet  ->fill(pT2,  weight);
	  if(momT>0.&&iPlus<3) {
	    _histZJet[iPlus]->fill(z,weight);
	    iPlus+=1;
	  }
	  else if(momT<0.&&iMinus<3) {
	    _histZJet[iMinus]->fill(z,weight);
	    iMinus+=1;
	  }
	}
      }
      _histPtTIn ->fill(pTSqIn /numParticles , weight);
      _histPtTOut->fill(pTSqOut/numParticles , weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_histSphericity);
      normalize(_histThrust);
      normalize(_histThrust2Jet);
      normalize(_histAplanarity);
      scale(_histZ        ,1./_wSum);
      scale(_histZScale   , sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      scale(_histXFeyn    ,1./_wSum/M_PI);
      scale(_histRap      ,1./_wSum);
      scale(_histZ2Jet    ,1./_wSum2);
      scale(_histXFeyn2Jet,1./_wSum2/M_PI);
      scale(_histRap2Jet  ,1./_wSum2);
      scale(_histPtT      ,1./_wSum);
      scale(_histPtT2Jet  ,1./_wSum2);
      scale(_histPtTIn    ,1./_wSum);
      scale(_histPtTOut   ,1./_wSum);
      for(unsigned int i=0;i<3;++i)
	scale(_histZJet[i]   ,0.5/_wSum2);
    }

    //@}


    /// @name Histograms
    //@{      
    Histo1DPtr _histSphericity, _histThrust, _histThrust2Jet, _histAplanarity,
      _histZ, _histZ2Jet, _histZScale, _histXFeyn, _histXFeyn2Jet, _histRap,
      _histRap2Jet, _histPtT, _histPtT2Jet, _histPtTIn, _histPtTOut ,_histZJet[3];
    double _wSum,_wSum2;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(HRS_1985_I201482);


}
