// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief pi, K and proton spectra at 12 and 30 GeV
  class TASSO_1980_I153656 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TASSO_1980_I153656);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");

      // Book histograms
      _iHist=-1;
      if(fuzzyEquals(sqrtS()/GeV, 12., 1e-3)) {
	_iHist = 0;
      }
      else if (fuzzyEquals(sqrtS()/GeV, 30., 1e-3)) {
	_iHist = 1;
      }
      else
	MSG_ERROR("Beam energy not supported!");

      _h_p_pi = bookHisto1D(3*_iHist+2,1,1);
      _h_x_pi = bookHisto1D(3*_iHist+2,1,2);
      _h_p_K  = bookHisto1D(3*_iHist+3,1,1);
      _h_x_K  = bookHisto1D(3*_iHist+3,1,2);
      _h_p_p  = bookHisto1D(3*_iHist+4,1,1);
      _h_x_p  = bookHisto1D(3*_iHist+4,1,2);

      _n_pi   = std::make_shared<YODA::Histo1D>(YODA::Histo1D(refData(3*_iHist+ 8,1,1)));
      _d_pi   = std::make_shared<YODA::Histo1D>(YODA::Histo1D(refData(3*_iHist+ 8,1,1)));
      _n_K    = std::make_shared<YODA::Histo1D>(YODA::Histo1D(refData(3*_iHist+ 9,1,1)));
      _d_K    = std::make_shared<YODA::Histo1D>(YODA::Histo1D(refData(3*_iHist+ 9,1,1)));
      _n_p    = std::make_shared<YODA::Histo1D>(YODA::Histo1D(refData(3*_iHist+10,1,1)));
      _d_p    = std::make_shared<YODA::Histo1D>(YODA::Histo1D(refData(3*_iHist+10,1,1)));

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const ChargedFinalState& fs = apply<ChargedFinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // Get event weight for histo filling
      const double weight = event.weight();

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      
      foreach (const Particle& p, fs.particles()) {
	double modp = p.p3().mod();
	_d_pi->fill(modp,weight);
	_d_K->fill(modp,weight);
	_d_p ->fill(modp,weight);
	double beta = modp/p.E();
	double xE = p.E()/meanBeamMom;
	if(abs(p.pdgId())==211) {
	  _h_p_pi->fill(modp, weight);
	  _h_x_pi->fill(xE  , weight/beta);
	  _n_pi->fill(modp,weight);
	}
	else if(abs(p.pdgId())==321) {
	  _h_p_K->fill(modp, weight);
	  _h_x_K->fill(xE  , weight/beta);
	  _n_K->fill(modp,weight);
	}
	else if(abs(p.pdgId())==2212) {
	  _h_p_p->fill(modp, weight);
	  _h_x_p->fill(xE  , weight/beta);
	  _n_p ->fill(modp,weight);
	}

      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_p_pi , crossSection()/nanobarn/sumOfWeights());
      scale(_h_x_pi , sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      scale(_h_p_K  , crossSection()/nanobarn/sumOfWeights());
      scale(_h_x_K  , sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      scale(_h_p_p  , crossSection()/nanobarn/sumOfWeights());
      scale(_h_x_p  , sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      
      divide(_n_pi,_d_pi, bookScatter2D(3*_iHist+ 8,1,1));
      divide(_n_K ,_d_K , bookScatter2D(3*_iHist+ 9,1,1));
      divide(_n_p ,_d_p , bookScatter2D(3*_iHist+10,1,1));
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_p_pi, _h_x_pi, _h_p_K, _h_x_K, _h_p_p, _h_x_p;
    Histo1DPtr _n_pi,_d_pi,_n_K,_d_K,_n_p,_d_p;
    int _iHist;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TASSO_1980_I153656);


}
