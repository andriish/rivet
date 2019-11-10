// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class TPC_1988_I262143 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TPC_1988_I262143);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");

      // Book histograms
      _h_z_pi   = bookHisto1D(1, 1, 1);
      _h_z_K    = bookHisto1D(1, 1, 2);
      _h_z_p    = bookHisto1D(1, 1, 3);
      _h_z_all  = bookHisto1D(1, 1, 4);
      
      _h_z2_pi   = bookHisto1D(5, 1, 1);
      _h_z2_K    = bookHisto1D(5, 1, 2);
      _h_z2_p    = bookHisto1D(5, 1, 3);
      
      _n_pi = std::make_shared<YODA::Histo1D>(Histo1D(refData(6,1,1)));
      _n_K  = std::make_shared<YODA::Histo1D>(Histo1D(refData(6,1,2)));
      _n_p  = std::make_shared<YODA::Histo1D>(Histo1D(refData(6,1,3)));
      _d_pi = std::make_shared<YODA::Histo1D>(Histo1D(refData(6,1,1)));
      _d_K  = std::make_shared<YODA::Histo1D>(Histo1D(refData(6,1,2)));
      _d_p  = std::make_shared<YODA::Histo1D>(Histo1D(refData(6,1,3)));
      _n2_K = std::make_shared<YODA::Histo1D>(Histo1D(refData(7,1,1)));
      _n2_p = std::make_shared<YODA::Histo1D>(Histo1D(refData(7,1,2)));
      _n3_p = std::make_shared<YODA::Histo1D>(Histo1D(refData(7,1,3)));
      _d2_K = std::make_shared<YODA::Histo1D>(Histo1D(refData(7,1,1)));
      _d2_p = std::make_shared<YODA::Histo1D>(Histo1D(refData(7,1,2)));
      _d3_p = std::make_shared<YODA::Histo1D>(Histo1D(refData(7,1,3)));
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

      // Get event weight for histo filling
      const double weight = event.weight();

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      foreach (const Particle& p, fs.particles()) {
	double xP = p.p3().mod()/meanBeamMom;
	_h_z_all->fill(xP, weight);
	_d_pi->fill(xP, weight);
	_d_K ->fill(xP, weight);
	_d_p ->fill(xP, weight);
	int id = abs(p.pdgId());
	if(id==211) {
	  _h_z_pi->fill(xP, weight);
	  _h_z2_pi->fill(xP, xP*weight);
	  _n_pi ->fill(xP, 100.*weight);
	  _d2_K->fill(xP, weight);
	  _d2_p->fill(xP, weight);
	  _d3_p->fill(xP, weight);
	}
	else if(id==321) {
	  _h_z_K ->fill(xP, weight);
	  _h_z2_K ->fill(xP, xP*weight);
	  _n_K ->fill(xP, 100.*weight);
	  _n2_K->fill(xP, weight);
	  _d3_p->fill(xP, weight);
	}
	else if(id==2212) {
	  _h_z_p ->fill(xP, weight);
	  _h_z2_p ->fill(xP, xP*weight);
	  _n_p  ->fill(xP, 100.*weight);
	  _n2_p->fill(xP, weight);
	  _n3_p->fill(xP, weight);
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_z_all,1./sumOfWeights());
      scale(_h_z_pi ,1./sumOfWeights());
      scale(_h_z_K  ,1./sumOfWeights());
      scale(_h_z_p  ,1./sumOfWeights());
      scale(_h_z2_pi ,1./sumOfWeights());
      scale(_h_z2_K  ,1./sumOfWeights());
      scale(_h_z2_p  ,1./sumOfWeights());
      divide(_n_pi,_d_pi, bookScatter2D(6, 1, 1));
      divide(_n_K ,_d_K , bookScatter2D(6, 1, 2));
      divide(_n_p ,_d_p , bookScatter2D(6, 1, 3));
      divide(_n2_K,_d2_K, bookScatter2D(7, 1, 1));
      divide(_n2_p,_d2_p, bookScatter2D(7, 1, 2));
      divide(_n3_p,_d3_p, bookScatter2D(7, 1, 3));
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_z_all,_h_z_pi,_h_z_K,_h_z_p;
    Histo1DPtr _h_z2_pi,_h_z2_K,_h_z2_p;
    Histo1DPtr _n_pi,_n_K,_n_p,_d_pi,_d_K,_d_p;
    Histo1DPtr _n2_K,_n2_p,_n3_p,_d2_K,_d2_p,_d3_p;
    //@}2


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TPC_1988_I262143);


}
