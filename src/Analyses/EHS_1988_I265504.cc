// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Event.hh"

namespace Rivet {


  class EHS_1988_I265504 : public Analysis {
  public:

    /// Constructor
    EHS_1988_I265504()
      : Analysis("EHS_1988_I265504")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// @todo Initialise and register projections here
      addProjection(ChargedFinalState(), "CFS");
      addProjection(Beam(),"Beam");

      /// @todo Book histograms here, e.g.:
      switch ( beamIds().first ) {
      case PID::PIPLUS:
	_h_cpos_xF = bookHisto1D(1, 1, 1);
	_h_cpos_eta = bookHisto1D(3, 1, 1);
	_h_cpos_pT2 = bookHisto1D(5, 1, 1);
	_h_cneg_xF = bookHisto1D(2, 1, 1);
	_h_cneg_eta = bookHisto1D(4, 1, 1);
	_h_cneg_pT2 = bookHisto1D(6, 1, 1);
	break;
	
      case PID::KPLUS:
	_h_cpos_xF = bookHisto1D(1, 1, 2);
	_h_cpos_eta = bookHisto1D(3, 1, 2);
	_h_cpos_pT2 = bookHisto1D(5, 1, 2);	
	_h_cneg_xF = bookHisto1D(2, 1, 2);
	_h_cneg_eta = bookHisto1D(4, 1, 2);
	_h_cneg_pT2 = bookHisto1D(6, 1, 2);
	break;
	
      case PID::PROTON:
	_h_cpos_xF = bookHisto1D(1, 1, 3);
	_h_cpos_eta = bookHisto1D(3, 1, 3);
	_h_cpos_pT2 = bookHisto1D(5, 1, 3);
	_h_cneg_xF = bookHisto1D(2, 1, 3);
       	_h_cneg_eta = bookHisto1D(4, 1, 3);
	_h_cneg_pT2 = bookHisto1D(6, 1, 3);	
	break;	
      }

      // calc. boost from lab to cm
      _beamboost = cmsTransform( beams() );
      MSG_INFO("boost vector: " << _beamboost );

      // transform beam into CMS frame
      Particle _beam_cm = beams().first;
      _beam_cm.transformBy(_beamboost);
      // beam momentum in CM frame defines Feynman-x
      _pz_max = _beam_cm.pz();

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
     
      const ChargedFinalState& fs = applyProjection<ChargedFinalState>(event, "CFS");
      Particles fs_particles = fs.particles();
      for (Particle& p: fs_particles ) {
	// only interested in pi- or positively charged
	if ( !( p.charge() > 0 || p.pid() == PID::PIMINUS ) ) continue;
	// slow proton cut: reject lab momenta < 1.2GeV
	if ( p.pid() == PID::PROTON && p.p()/GeV < 1.2  ) continue;
	// transform to cm frame
	p.transformBy(_beamboost);	
	const double xF = p.pz()/_pz_max;
	
	if( p.charge() > 0 ){	
	  _h_cpos_xF->fill( xF, weight);	  
	  _h_cpos_pT2->fill( p.pT2(), weight);
	  _h_cpos_eta->fill( p.eta(), weight);
	} else {
	  if ( p.pid() == PID::PIMINUS ) {
	    _h_cneg_xF->fill( xF, weight);
	    _h_cneg_pT2->fill( p.pT2(), weight);
	    _h_cneg_eta->fill( p.eta(), weight);
	  } else {
	    continue;
	  }	  
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // , _h_cneg_eta,  _h_cneg_pT2
      scale( {_h_cpos_xF, _h_cpos_pT2,_h_cpos_eta, _h_cneg_xF, _h_cneg_eta, _h_cneg_pT2 } , crossSection()/millibarn/sumOfWeights() );

    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here


    /// @name Histograms
    //@{
    LorentzTransform _beamboost;
    double _pz_max;
    Histo1DPtr _h_cpos_xF;
    Histo1DPtr _h_cpos_eta;
    Histo1DPtr _h_cpos_pT2;

    Histo1DPtr _h_cneg_xF;
    Histo1DPtr _h_cneg_eta;
    Histo1DPtr _h_cneg_pT2;

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(EHS_1988_I265504);


}
