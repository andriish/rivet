// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"

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
      const FinalState fs;
      declare(fs, "FS");
      const Thrust thrust(fs);
      declare(thrust, "Thrust");
      declare(Sphericity(fs), "Sphericity");
      declare(ParisiTensor(fs), "Parisi");
      declare(Hemispheres(thrust), "Hemispheres");

      // find the beam energy
      unsigned int offset = 0;
      if      (fuzzyEquals(sqrtS()/GeV, 45, 1E-3)) offset = 1;
      else if (fuzzyEquals(sqrtS()/GeV, 66, 1E-3)) offset = 2;
      else if (fuzzyEquals(sqrtS()/GeV, 76, 1E-3)) offset = 3;
      else    MSG_ERROR("Beam energy not supported!");
      // Book the histograms
      book(_h_thrust,  1, 1, offset);
      book(_h_major,  2, 1, offset);
      book(_h_minor,  3, 1, offset);
      book(_h_sphericity,  4, 1, offset);
      book(_h_planarity,  5, 1, offset);
      book(_h_oblateness,  6, 1, offset);
      book(_h_heavy_jet_mass,  7, 1, offset);
      book(_h_light_jet_mass,  9, 1, offset);
      book(_h_diff_jet_mass, 10, 1, offset);
      book(_h_total_jet_mass, 11, 1, offset);
      book(_h_heavy_jet_mass_E,  8, 1, offset);
      book(_h_total_jet_mass_E, 12, 1, offset);
      book(_h_wide_broading, 13, 1, offset);
      book(_h_narrow_broading, 14, 1, offset);
      book(_h_total_broading, 15, 1, offset);
      book(_h_diff_broading, 16, 1, offset);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      // thrust related observables
      _h_thrust    ->fill(1.-thrust.thrust()  );
      _h_major     ->fill(thrust.thrustMajor());
      _h_minor     ->fill(thrust.thrustMinor());
      _h_oblateness->fill(thrust.oblateness() );

      // sphericity related
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      _h_sphericity->fill(sphericity.sphericity());
      _h_planarity ->fill(sphericity.planarity() );
      // hemisphere related
      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");
      // standard jet masses
      _h_heavy_jet_mass->fill(hemi.scaledM2high());
      _h_light_jet_mass->fill(hemi.scaledM2low() );
      _h_diff_jet_mass ->fill(hemi.scaledM2diff());
      _h_total_jet_mass->fill(hemi.scaledM2low()+hemi.scaledM2high());
      // jet broadening
      _h_wide_broading  ->fill(hemi.Bmax() );
      _h_narrow_broading->fill(hemi.Bmin() );
      _h_total_broading ->fill(hemi.Bsum() );
      _h_diff_broading  ->fill(hemi.Bdiff());
      // E scheme jet masses
      Vector3 axis = thrust.thrustAxis();
      FourMomentum p4With, p4Against;
      double Evis(0);
      for (const Particle& p : apply<FinalState>(event, "FS").particles()) {
	Vector3 p3 = p.momentum().vector3().unitVec();
	const double   E = p.momentum().E();
	Evis += E;
	p3 = E*p3;
	const double p3Para = dot(p3, axis);
	FourMomentum p4(E,p3.x(),p3.y(),p3.z());
	if (p3Para > 0)      p4With    += p4;
	else if (p3Para < 0) p4Against += p4;
	else {
	  MSG_WARNING("Particle split between hemispheres");
	  p4With    += 0.5 * p4;
	  p4Against += 0.5 * p4;
	}
      }
      const double mass2With    = p4With.mass2()/sqr(Evis);
      const double mass2Against = p4Against.mass2()/sqr(Evis);
      // fill the histograms
      _h_heavy_jet_mass_E->fill(max(mass2With,mass2Against));
      _h_total_jet_mass_E->fill(mass2With+mass2Against);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_thrust          );
      normalize(_h_major           );
      normalize(_h_minor           );
      normalize(_h_sphericity      );
      normalize(_h_planarity       );
      normalize(_h_oblateness      );
      normalize(_h_heavy_jet_mass  );
      normalize(_h_light_jet_mass  );
      normalize(_h_diff_jet_mass   );
      normalize(_h_total_jet_mass  );
      normalize(_h_heavy_jet_mass_E);
      normalize(_h_total_jet_mass_E);
      normalize(_h_wide_broading   );
      normalize(_h_narrow_broading );
      normalize(_h_total_broading  );
      normalize(_h_diff_broading   );

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_thrust,_h_major,_h_minor;
    Histo1DPtr _h_sphericity,_h_planarity,_h_oblateness;
    Histo1DPtr _h_heavy_jet_mass,_h_light_jet_mass,_h_diff_jet_mass,_h_total_jet_mass;
    Histo1DPtr _h_heavy_jet_mass_E,_h_total_jet_mass_E;
    Histo1DPtr _h_wide_broading,_h_narrow_broading,_h_total_broading,_h_diff_broading;
    //@}
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(DELPHI_2003_I620250);


}
