// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief K0 and Lambda spectra at 14, 22 and 34 GeC
  class TASSO_1985_I205119 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TASSO_1985_I205119);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      // Book histograms
      if(fuzzyEquals(sqrtS()/GeV, 14., 1e-3)) {
	_h_kaon_x   = bookHisto1D( 1,1,1);
	_h_lambda_x = bookHisto1D( 4,1,1);
	_h_kaon_p   = bookHisto1D( 7,1,1);
	_h_lambda_p = bookHisto1D(10,1,1);
      }
      else if (fuzzyEquals(sqrtS()/GeV, 22., 1e-3)) {
	_h_kaon_x   = bookHisto1D( 2,1,1);
	_h_lambda_x = bookHisto1D( 5,1,1);
	_h_kaon_p   = bookHisto1D( 8,1,1);
	_h_lambda_p = bookHisto1D(11,1,1);
      }
      else if (fuzzyEquals(sqrtS()/GeV, 34., 1e-3)) {
	_h_kaon_x   = bookHisto1D( 3,1,1);
	_h_lambda_x = bookHisto1D( 6,1,1);
	_h_kaon_p   = bookHisto1D( 9,1,1);
	_h_lambda_p = bookHisto1D(12,1,1);
      }
      else
	MSG_ERROR("Beam energy not supported!");
	
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

      foreach (const Particle& p, apply<UnstableParticles>(event, "UFS").
	       particles(Cuts::abspid==PID::LAMBDA or Cuts::pid==130 or Cuts::pid==310)) {
	double xE = p.E()/meanBeamMom;
	double modp = p.p3().mod();
	double beta = modp/p.E();
	if(abs(p.pdgId())==PID::LAMBDA) {
	  _h_lambda_x->fill(xE,weight/beta);
	  _h_lambda_p->fill(modp,weight);
	}
	else {
	  _h_kaon_x->fill(xE,weight/beta);	 
	  _h_kaon_p->fill(modp,weight);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_kaon_x  , sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      scale(_h_lambda_x, sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      scale(_h_kaon_p  , crossSection()/nanobarn/sumOfWeights());
      scale(_h_lambda_p, crossSection()/nanobarn/sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_kaon_x, _h_lambda_x, _h_kaon_p, _h_lambda_p;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TASSO_1985_I205119);


}
