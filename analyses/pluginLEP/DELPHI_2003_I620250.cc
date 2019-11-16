// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/ParisiTensor.hh"

namespace Rivet {


  /// @brief DELPHI event shapes below the Z pole
  class DELPHI_2003_I620250 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(DELPHI_2003_I620250);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      const FinalState fs;
      declare(fs, "FS");
      const Thrust thrust(fs);
      declare(thrust, "Thrust");
      declare(Sphericity(fs), "Sphericity");
      declare(ParisiTensor(fs), "Parisi");
      declare(Hemispheres(thrust), "Hemispheres");

      // find the beam energy
      unsigned int offset = 0;
      int offset2 = -1;
      if      (fuzzyEquals(sqrtS()/GeV, 45, 1E-3)) offset = 1;
      else if (fuzzyEquals(sqrtS()/GeV, 66, 1E-3)) offset = 2;
      else if (fuzzyEquals(sqrtS()/GeV, 76, 1E-3)) offset = 3;
      else if (fuzzyEquals(sqrtS()/GeV, 183  , 1E-3)) {
	offset2= 0;			   
	offset = 1;			   
      }					   
      else if (fuzzyEquals(sqrtS()/GeV, 189  , 1E-3)) {
	offset2= 0;			   
	offset = 2;			   
      }					   
      else if (fuzzyEquals(sqrtS()/GeV, 192  , 1E-3)) {
	offset2= 0;			   
	offset = 3;			   
      }					   
      else if (fuzzyEquals(sqrtS()/GeV, 196  , 1E-3)) {
	offset2= 0;			   
	offset = 4;			   
      }					   
      else if (fuzzyEquals(sqrtS()/GeV, 200  , 1E-3)) {
	offset2= 1;			   
	offset = 1;			   
      }					   
      else if (fuzzyEquals(sqrtS()/GeV, 202  , 1E-3)) {
	offset2= 1;			   
	offset = 2;			   
      }					   
      else if (fuzzyEquals(sqrtS()/GeV, 205  , 1E-3)) {
	offset2= 1;			   
	offset = 3;			   
      }					   
      else if (fuzzyEquals(sqrtS()/GeV, 207  , 1E-3)) {
	offset2= 1;
	offset = 4;
      }
      else    MSG_ERROR("Beam energy not supported!");
      // Book the histograms
      if(offset2<0) {
	_h_thrust           = bookHisto1D( 1, 1, offset);
	_h_major            = bookHisto1D( 2, 1, offset);
	_h_minor            = bookHisto1D( 3, 1, offset);
	_h_sphericity       = bookHisto1D( 4, 1, offset);
	_h_planarity        = bookHisto1D( 5, 1, offset);
	_h_oblateness       = bookHisto1D( 6, 1, offset);
	_h_heavy_jet_mass   = bookHisto1D( 7, 1, offset);
	_h_light_jet_mass   = bookHisto1D( 9, 1, offset);
	_h_diff_jet_mass    = bookHisto1D(10, 1, offset);
	_h_total_jet_mass   = bookHisto1D(11, 1, offset);
	_h_heavy_jet_mass_E = bookHisto1D( 8, 1, offset);
	_h_total_jet_mass_E = bookHisto1D(12, 1, offset);
	_h_wide_broading    = bookHisto1D(13, 1, offset);
	_h_narrow_broading  = bookHisto1D(14, 1, offset);
	_h_total_broading   = bookHisto1D(15, 1, offset);
	_h_diff_broading    = bookHisto1D(16, 1, offset);
	_h_CParam           = bookHisto1D(17, 1, offset);
      }
      else {
	_h_rap              = bookHisto1D(30+offset2, 1, offset);
	_h_xi               = bookHisto1D(32+offset2, 1, offset);
	_h_pTIn             = bookHisto1D(34+offset2, 1, offset);
	_h_pTOut            = bookHisto1D(36+offset2, 1, offset);
	_h_thrust           = bookHisto1D(38+offset2, 1, offset);
	_h_major            = bookHisto1D(40+offset2, 1, offset);
	_h_minor            = bookHisto1D(42+offset2, 1, offset);
	_h_oblateness       = bookHisto1D(44+offset2, 1, offset);
	_h_wide_broading    = bookHisto1D(46+offset2, 1, offset);
	_h_total_broading   = bookHisto1D(48+offset2, 1, offset);
	_h_diff_broading    = bookHisto1D(50+offset2, 1, offset);
	_h_CParam           = bookHisto1D(52+offset2, 1, offset);
	_h_DParam           = bookHisto1D(54+offset2, 1, offset);
	_h_heavy_jet_mass   = bookHisto1D(56+offset2, 1, offset);
	_h_heavy_jet_mass_P = bookHisto1D(58+offset2, 1, offset);
	_h_heavy_jet_mass_E = bookHisto1D(60+offset2, 1, offset);
	_h_light_jet_mass   = bookHisto1D(62+offset2, 1, offset);
	_h_diff_jet_mass    = bookHisto1D(64+offset2, 1, offset);
	_h_sphericity       = bookHisto1D(66+offset2, 1, offset);
	_h_planarity        = bookHisto1D(68+offset2, 1, offset);
	_h_aplanarity       = bookHisto1D(70+offset2, 1, offset);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      // thrust related observables
      _h_thrust    ->fill(1.-thrust.thrust()  ,weight);
      _h_major     ->fill(thrust.thrustMajor(),weight);
      _h_minor     ->fill(thrust.thrustMinor(),weight);
      _h_oblateness->fill(thrust.oblateness() ,weight);

      // sphericity related
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      _h_sphericity->fill(sphericity.sphericity(),weight);
      _h_planarity ->fill(sphericity.planarity() ,weight);
      if(_h_aplanarity) _h_aplanarity->fill(sphericity.aplanarity() ,weight);
      // hemisphere related
      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");
      // standard jet masses
      _h_heavy_jet_mass->fill(hemi.scaledM2high(),weight);
      _h_light_jet_mass->fill(hemi.scaledM2low() ,weight);
      _h_diff_jet_mass ->fill(hemi.scaledM2diff(),weight);
      if(_h_total_jet_mass) _h_total_jet_mass->fill(hemi.scaledM2low()+hemi.scaledM2high(),weight);
      // jet broadening
      _h_wide_broading  ->fill(hemi.Bmax() ,weight);
      if(_h_narrow_broading) _h_narrow_broading->fill(hemi.Bmin() ,weight);
      _h_total_broading ->fill(hemi.Bsum() ,weight);
      _h_diff_broading  ->fill(hemi.Bdiff(),weight);
      // E and p scheme jet masses
      Vector3 axis = thrust.thrustAxis();
      FourMomentum p4WithE, p4AgainstE;
      FourMomentum p4WithP, p4AgainstP;
      double Evis(0);
      foreach(const Particle& p, apply<FinalState>(event, "FS").particles()) {
	Vector3 p3 = p.momentum().vector3().unitVec();
	const double   E = p.momentum().E();
	Evis += E;
	p3 = E*p3;
	const double p3Para = dot(p3, axis);
	FourMomentum p4E(E,p3.x(),p3.y(),p3.z());
	FourMomentum p4P(p.p3().mod(),p.p3().x(),p.p3().y(),p.p3().z());
	if (p3Para > 0)      {
	  p4WithE    += p4E;
	  p4WithP    += p4P;
	}
	else if (p3Para < 0) {
	  p4AgainstE += p4E;
	  p4AgainstP += p4P;
	}
	else {
	  MSG_WARNING("Particle split between hemispheres");
	  p4WithE    += 0.5 * p4E;
	  p4AgainstE += 0.5 * p4E;
	  p4WithP    += 0.5 * p4P;
	  p4AgainstP += 0.5 * p4P;
	}
      }
      // E scheme
      const double mass2With_E    = p4WithE.mass2()/sqr(Evis);
      const double mass2Against_E = p4AgainstE.mass2()/sqr(Evis);
      // fill the histograms
      _h_heavy_jet_mass_E->fill(max(mass2With_E,mass2Against_E),weight);
      if(_h_total_jet_mass_E) _h_total_jet_mass_E->fill(mass2With_E+mass2Against_E,weight);
      // pscheme
      const double mass2With_P    = p4WithP.mass2()/sqr(Evis);
      const double mass2Against_P = p4AgainstP.mass2()/sqr(Evis);
      // fill the histograms
      if(_h_heavy_jet_mass_P) _h_heavy_jet_mass_P->fill(max(mass2With_P,mass2Against_P),weight);
      
      MSG_DEBUG("Calculating Parisi params");
      const ParisiTensor& parisi = apply<ParisiTensor>(event, "Parisi");
      _h_CParam->fill(parisi.C(), weight);
      if(_h_DParam) _h_DParam->fill(parisi.D(), weight);

      // single particle distributions
      const FinalState& fs = apply<FinalState>(event, "FS");
      if(_h_xi) {
	foreach (const Particle& p, fs.particles()) {
	  if( ! PID::isCharged(p.pdgId())) continue;
	  // Get momentum and energy of each particle.
	  const Vector3 mom3 = p.p3();
	  const double energy = p.E();
	  
	  // Scaled momenta.
	  const double mom = mom3.mod();
	  const double scaledMom = mom/meanBeamMom;
	  const double logInvScaledMom = -std::log(scaledMom);
	  _h_xi->fill(logInvScaledMom, weight);
	  
	  // Get momenta components w.r.t. thrust and sphericity.
	  const double momT = dot(thrust.thrustAxis(), mom3);
	  const double pTinT = dot(mom3, thrust.thrustMajorAxis());
	  const double pToutT = dot(mom3, thrust.thrustMinorAxis());
	  _h_pTIn ->fill(fabs(pTinT/GeV), weight);
	  _h_pTOut->fill(fabs(pToutT/GeV), weight);
	  
	  // Calculate rapidities w.r.t. thrust and sphericity.
	  const double rapidityT = 0.5 * std::log((energy + momT) / (energy - momT));
	  _h_rap->fill(fabs(rapidityT), weight);
	  MSG_TRACE(fabs(rapidityT) << " " << scaledMom/GeV);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_thrust          );
      normalize(_h_major           );
      normalize(_h_minor           );
      normalize(_h_sphericity      );
      normalize(_h_planarity       );
      if(_h_aplanarity) normalize(_h_aplanarity       );
      normalize(_h_oblateness      );
      normalize(_h_heavy_jet_mass  );
      normalize(_h_light_jet_mass  );
      normalize(_h_diff_jet_mass   );
      if(_h_total_jet_mass) normalize(_h_total_jet_mass  );
      normalize(_h_heavy_jet_mass_E);
      if(_h_total_jet_mass_E) normalize(_h_total_jet_mass_E);
      if(_h_heavy_jet_mass_P) normalize(_h_heavy_jet_mass_P);
      normalize(_h_wide_broading   );
      if(_h_narrow_broading) normalize(_h_narrow_broading );
      normalize(_h_total_broading  );
      normalize(_h_diff_broading   );
      normalize(_h_CParam   );
      if(_h_DParam) normalize(_h_DParam   );
      if(_h_xi) {
	scale(_h_xi   ,1./sumOfWeights());
	scale(_h_pTIn ,1./sumOfWeights());
	scale(_h_pTOut,1./sumOfWeights());
	scale(_h_rap  ,1./sumOfWeights());
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_thrust,_h_major,_h_minor;
    Histo1DPtr _h_sphericity,_h_planarity,_h_aplanarity,_h_oblateness;
    Histo1DPtr _h_heavy_jet_mass,_h_light_jet_mass,_h_diff_jet_mass,_h_total_jet_mass;
    Histo1DPtr _h_heavy_jet_mass_E,_h_total_jet_mass_E;
    Histo1DPtr _h_heavy_jet_mass_P;
    Histo1DPtr _h_wide_broading,_h_narrow_broading,_h_total_broading,_h_diff_broading;
    Histo1DPtr _h_CParam,_h_DParam;
    Histo1DPtr _h_xi, _h_pTIn, _h_pTOut,_h_rap;
    //@}
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(DELPHI_2003_I620250);


}
