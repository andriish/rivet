// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InitialQuarks.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"

namespace Rivet {


  /// Jet rates and event shapes at LEP I+II
  class L3_2004_I652683 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(L3_2004_I652683);

    /// Book histograms and initialise projections before the run
    void init() {

      // Projections to use
      const FinalState FS;
      declare(FS, "FS");
      declare(Beam(), "beams");
      const ChargedFinalState CFS;
      declare(CFS, "CFS");
      const Thrust thrust(FS);
      declare(thrust, "thrust");
      declare(ParisiTensor(FS), "Parisi");
      declare(Hemispheres(thrust), "Hemispheres");
      declare(InitialQuarks(), "initialquarks");

      // Book the histograms
      if(fuzzyEquals(sqrtS()/GeV, 91.2, 1e-3)) {
	// z pole
	_h_Thrust_udsc             = bookHisto1D(47, 1, 1);
	_h_Thrust_bottom           = bookHisto1D(47, 1, 2);
	_h_heavyJetmass_udsc       = bookHisto1D(48, 1, 1);
	_h_heavyJetmass_bottom     = bookHisto1D(48, 1, 2);
	_h_totalJetbroad_udsc      = bookHisto1D(49, 1, 1);
	_h_totalJetbroad_bottom    = bookHisto1D(49, 1, 2);
	_h_wideJetbroad_udsc       = bookHisto1D(50, 1, 1);
	_h_wideJetbroad_bottom     = bookHisto1D(50, 1, 2);
	_h_Cparameter_udsc         = bookHisto1D(51, 1, 1);
	_h_Cparameter_bottom       = bookHisto1D(51, 1, 2);
	_h_Dparameter_udsc         = bookHisto1D(52, 1, 1);
	_h_Dparameter_bottom       = bookHisto1D(52, 1, 2);
	_h_Ncharged                = bookHisto1D("/TMP/NCHARGED"     , 28, 1, 57);
	_h_Ncharged_udsc           = bookHisto1D("/TMP/NCHARGED_UDSC", 28, 1, 57);
	_h_Ncharged_bottom         = bookHisto1D("/TMP/NCHARGED_B"   , 27, 3, 57);
	_h_scaledMomentum          = bookHisto1D(65, 1, 1);
	_h_scaledMomentum_udsc     = bookHisto1D(65, 1, 2);
	_h_scaledMomentum_bottom   = bookHisto1D(65, 1, 3);
      }
      else if(sqrtS()/GeV<90) {
	int i1(-1),i2(-1);
	if(fuzzyEquals(sqrtS()/GeV, 41.4, 1e-2)) {
	  i1=0;
	  i2=1;
	}
	else if(fuzzyEquals(sqrtS()/GeV, 55.3, 1e-2)) {
	  i1=0;
	  i2=2;
	}
	else if(fuzzyEquals(sqrtS()/GeV, 65.4, 1e-2)) {
	  i1=0;
	  i2=3;
	}
	else if(fuzzyEquals(sqrtS()/GeV, 75.7, 1e-2)) {
	  i1=1;
	  i2=1;
	}
	else if(fuzzyEquals(sqrtS()/GeV, 82.3, 1e-2)) {
	  i1=1;
	  i2=2;
	}
	else if(fuzzyEquals(sqrtS()/GeV, 85.1, 1e-2)) {
	  i1=1;
	  i2=3;
	}
	else
	  MSG_ERROR("Beam energy not supported!");
	_h_thrust = bookHisto1D(21+i1,1,i2);
	_h_rho    = bookHisto1D(26+i1,1,i2);
	_h_B_T    = bookHisto1D(31+i1,1,i2);
	_h_B_W    = bookHisto1D(36+i1,1,i2);
      }
      else if(sqrtS()/GeV>120) {
	int i1(-1),i2(-1);
	if(fuzzyEquals(sqrtS()/GeV, 130.1, 1e-2)) {
	  i1=0;
	  i2=1;
	}
	else if(fuzzyEquals(sqrtS()/GeV, 136.1, 1e-2)) {
	  i1=0;
	  i2=2;
	}
	else if(fuzzyEquals(sqrtS()/GeV, 161.3, 1e-2)) {
	  i1=0;
	  i2=3;
	}
	else if(fuzzyEquals(sqrtS()/GeV, 172.3, 1e-2)) {
	  i1=1;
	  i2=1;
	}
	else if(fuzzyEquals(sqrtS()/GeV, 182.8, 1e-2)) {
	  i1=1;
	  i2=2;
	}
	else if(fuzzyEquals(sqrtS()/GeV, 188.6, 1e-2)) {
	  i1=1;
	  i2=3;
	}
	else if(fuzzyEquals(sqrtS()/GeV, 194.4, 1e-2)) {
	  i1=2;
	  i2=1;
	}
	else if(fuzzyEquals(sqrtS()/GeV, 200.2, 1e-2)) {
	  i1=2;
	  i2=2;
	}
	else if(fuzzyEquals(sqrtS()/GeV, 206.2, 1e-2)) {
	  i1=2;
	  i2=3;
	}
	else
	  MSG_ERROR("Beam energy not supported!");
	_h_thrust = bookHisto1D(23+i1,1,i2);
	_h_rho    = bookHisto1D(28+i1,1,i2);
	_h_B_T    = bookHisto1D(33+i1,1,i2);
	_h_B_W    = bookHisto1D(38+i1,1,i2);
	_h_C      = bookHisto1D(41+i1,1,i2);
	_h_D      = bookHisto1D(44+i1,1,i2);
	_h_N      = bookHisto1D("/TMP/NCHARGED", 22, 9, 53);
	_h_xi     = bookHisto1D(66+i1,1,i2);
	// todo add the jets
	// int i3 = 3*i1+i2;
	// _h_y_2_JADE   = bookHisto1D(   i3,1,1);
	// _h_y_3_JADE   = bookHisto1D(   i3,1,2);
	// _h_y_4_JADE   = bookHisto1D(   i3,1,3);
	// _h_y_5_JADE   = bookHisto1D(   i3,1,4);
	// _h_y_2_Durham = bookHisto1D( 9+i3,1,1);
	// _h_y_3_Durham = bookHisto1D( 9+i3,1,2);
	// _h_y_4_Durham = bookHisto1D( 9+i3,1,3);
	// _h_y_5_Durham = bookHisto1D( 9+i3,1,4);
	// if(i3==8||i3==9) {
	//   _h_y_2_Cambridge = bookHisto1D(10+i3,1,1);
	//   _h_y_3_Cambridge = bookHisto1D(10+i3,1,2);
	//   _h_y_4_Cambridge = bookHisto1D(10+i3,1,3);
	//   _h_y_5_Cambridge = bookHisto1D(10+i3,1,4);
	// }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beam average momentum
      const ParticlePair& beams = apply<Beam>(event, "beams").beams();
      const double beamMomentum = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;

      // InitialQuarks projection to have udsc events separated from b events
      /// @todo Yuck!!! Eliminate when possible...
      int iflav = 0;
      // only need the flavour at Z pole
      if(_h_Thrust_udsc) {
	int flavour = 0;
	const InitialQuarks& iqf = apply<InitialQuarks>(event, "initialquarks");
	Particles quarks;
	if ( iqf.particles().size() == 2 ) {
	  flavour = iqf.particles().front().abspid();
	  quarks  = iqf.particles();
	} else {
	  map<int, Particle> quarkmap;
	  for (const Particle& p : iqf.particles()) {
	    if (quarkmap.find(p.pid()) == quarkmap.end()) quarkmap[p.pid()] = p;
	    else if (quarkmap[p.pid()].E() < p.E()) quarkmap[p.pid()] = p;
	  }
	  double max_energy = 0.;
	  for (int i = 1; i <= 5; ++i) {
	    double energy = 0.;
	    if (quarkmap.find(i) != quarkmap.end())
	      energy += quarkmap[ i].E();
	    if (quarkmap.find(-i) != quarkmap.end())
	      energy += quarkmap[-i].E();
	    if (energy > max_energy)
	      flavour = i;
	  }
	  if (quarkmap.find(flavour) != quarkmap.end())
	    quarks.push_back(quarkmap[flavour]);
	  if (quarkmap.find(-flavour) != quarkmap.end())
	    quarks.push_back(quarkmap[-flavour]);
	}
	// Flavour label
	/// @todo Change to a bool?
	iflav = (flavour == PID::DQUARK || flavour == PID::UQUARK || flavour == PID::SQUARK || flavour == PID::CQUARK) ? 1 : (flavour == PID::BQUARK) ? 5 : 0;
      }
      // Update weight sums
      const double weight = event.weight();
      if (iflav == 1) {
        _sumW_udsc += weight;
      } else if (iflav == 5) {
        _sumW_b += weight;
      }
      _sumW_ch += weight;
      
      // Charged multiplicity
      const FinalState& cfs = applyProjection<FinalState>(event, "CFS");
      if(_h_Ncharged) _h_Ncharged->fill(cfs.size(), weight);
      if (iflav == 1) {
        _sumW_ch_udsc += weight;
        _h_Ncharged_udsc->fill(cfs.size(), weight);
      }
      else if (iflav == 5) {
        _sumW_ch_b += weight;
        _h_Ncharged_bottom->fill(cfs.size(), weight);
      }
      else if(_h_N) {
	_h_N->fill(cfs.size(),weight);
      }

      // Scaled momentum
      const Particles& chparticles = cfs.particlesByPt();
      for (const Particle& p : chparticles) {
        const Vector3 momentum3 = p.p3();
        const double mom = momentum3.mod();
        const double scaledMom = mom/beamMomentum;
        const double logScaledMom = std::log(scaledMom);
        if(_h_scaledMomentum) _h_scaledMomentum->fill(-logScaledMom, weight);
        if (iflav == 1) {
          _h_scaledMomentum_udsc->fill(-logScaledMom, weight);
        } else if (iflav == 5) {
          _h_scaledMomentum_bottom->fill(-logScaledMom, weight);
        }
	else if(_h_xi) {
	  _h_xi->fill(-logScaledMom, weight);
	}
      }
      
      // Thrust
      const Thrust& thrust = applyProjection<Thrust>(event, "thrust");
      if (iflav == 1) {
        _h_Thrust_udsc->fill(thrust.thrust(), weight);
      }
      else if (iflav == 5) {
        _h_Thrust_bottom->fill(thrust.thrust(), weight);
      }
      else if(_h_thrust) {
        _h_thrust->fill(1.-thrust.thrust(), weight);
      }

      // C and D Parisi parameters
      const ParisiTensor& parisi = applyProjection<ParisiTensor>(event, "Parisi");
      if (iflav == 1) {
        _h_Cparameter_udsc->fill(parisi.C(), weight);
        _h_Dparameter_udsc->fill(parisi.D(), weight);
      } else if (iflav == 5) {
        _h_Cparameter_bottom->fill(parisi.C(), weight);
        _h_Dparameter_bottom->fill(parisi.D(), weight);
      }
      else if(_h_C) {
	_h_C->fill(parisi.C(), weight);
	_h_D->fill(parisi.D(), weight);
      }

      // The hemisphere variables
      const Hemispheres& hemisphere = applyProjection<Hemispheres>(event, "Hemispheres");
      if (iflav == 1) {
        _h_heavyJetmass_udsc->fill(hemisphere.scaledM2high(), weight);
        _h_totalJetbroad_udsc->fill(hemisphere.Bsum(), weight);
        _h_wideJetbroad_udsc->fill(hemisphere.Bmax(), weight);
      }
      else if (iflav == 5) {
        _h_heavyJetmass_bottom->fill(hemisphere.scaledM2high(), weight);
        _h_totalJetbroad_bottom->fill(hemisphere.Bsum(), weight);
        _h_wideJetbroad_bottom->fill(hemisphere.Bmax(), weight);
      }
      else if (_h_rho) {
        _h_rho->fill(hemisphere.scaledM2high(), weight);
        _h_B_T->fill(hemisphere.Bsum(), weight);
        _h_B_W->fill(hemisphere.Bmax(), weight);
      }

    }
    
    Scatter2DPtr convertHisto(unsigned int ix,unsigned int iy, unsigned int iz, Histo1DPtr histo) {
      Scatter2D temphisto(refData(ix, iy, iz));
      Scatter2DPtr mult = bookScatter2D(ix, iy, iz);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
	const double x  = temphisto.point(b).x();
	pair<double,double> ex = temphisto.point(b).xErrs();
	double y    = histo->bins()[b].area();
	double yerr = histo->bins()[b].areaErr();
	mult->addPoint(x, y, ex, make_pair(yerr,yerr));
      }
      return mult;
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // Z pole plots
      if(_h_Thrust_udsc) {
	scale({_h_Thrust_udsc, _h_heavyJetmass_udsc, _h_totalJetbroad_udsc,
	      _h_wideJetbroad_udsc, _h_Cparameter_udsc, _h_Dparameter_udsc}, 1/_sumW_udsc);
	scale({_h_Thrust_bottom, _h_heavyJetmass_bottom, _h_totalJetbroad_bottom,
	      _h_wideJetbroad_bottom, _h_Cparameter_bottom, _h_Dparameter_bottom}, 1./_sumW_b);
	scale(_h_Ncharged, 1./_sumW_ch);
	scale(_h_Ncharged_udsc, 1./_sumW_ch_udsc);
	scale(_h_Ncharged_bottom, 1./_sumW_ch_b);
	convertHisto(59,1,1,_h_Ncharged       );
	convertHisto(59,1,2,_h_Ncharged_udsc  );
	convertHisto(59,1,3,_h_Ncharged_bottom);


	
	scale(_h_scaledMomentum, 1/_sumW_ch);
	scale(_h_scaledMomentum_udsc, 1/_sumW_ch_udsc);
	scale(_h_scaledMomentum_bottom, 1/_sumW_ch_b);
      }
      else {
	if(_h_thrust) normalize(_h_thrust);
	if(_h_rho) normalize(_h_rho);
	if(_h_B_T) normalize(_h_B_T);
	if(_h_B_W) normalize(_h_B_W);
	if(_h_C) normalize(_h_C);
	if(_h_D) normalize(_h_D);
	if(_h_N) normalize(_h_N);
	if(_h_xi) scale(_h_xi,1./sumOfWeights());
	
      
      Scatter2DPtr mult;
	if(_h_N) {
	  if(fuzzyEquals(sqrtS()/GeV, 130.1, 1e-2)) {
	    convertHisto(60, 1, 1, _h_N);
	  }
	  else if(fuzzyEquals(sqrtS()/GeV, 136.1, 1e-2)) {
	    convertHisto(60, 1, 2, _h_N);
	  }
	  else if(fuzzyEquals(sqrtS()/GeV, 161.3, 1e-2)) {
	    convertHisto(60, 1, 3, _h_N);
	  }
	  else if(fuzzyEquals(sqrtS()/GeV, 172.3, 1e-2)) {
	    convertHisto(61, 1, 1, _h_N);
	  }
	  else if(fuzzyEquals(sqrtS()/GeV, 182.8, 1e-2)) {
	    convertHisto(61, 1, 2, _h_N);
	  }
	  else if(fuzzyEquals(sqrtS()/GeV, 188.6, 1e-2)) {
	    convertHisto(61, 1, 3, _h_N);
	  }
	  else if(fuzzyEquals(sqrtS()/GeV, 194.4, 1e-2)) {
	    convertHisto(62, 1, 1, _h_N);
	  }
	  else if(fuzzyEquals(sqrtS()/GeV, 200.2, 1e-2)) {
	    convertHisto(62, 1, 2, _h_N);
	  }
	  else if(fuzzyEquals(sqrtS()/GeV, 206.2, 1e-2)) {
	    convertHisto(62, 1, 3, _h_N);
	  }
	}
	// todo add the jets
	// Histo1DPtr _h_y_2_JADE,_h_y_3_JADE,_h_y_4_JADE,_h_y_5_JADE;
	// Histo1DPtr _h_y_2_Durham,_h_y_3_Durham,_h_y_4_Durham,_h_y_5_Durham;
	// Histo1DPtr _h_y_2_Cambridge,_h_y_3_Cambridge,_h_y_4_Cambridge,_h_y_5_Cambridge;
      }
    }


    /// Weight counters
    double _sumW_udsc = 0, _sumW_b = 0, _sumW_ch = 0, _sumW_ch_udsc = 0, _sumW_ch_b = 0;

    /// @name Histograms
    //@{
    // at the Z pole
    Histo1DPtr _h_Thrust_udsc, _h_Thrust_bottom;
    Histo1DPtr _h_heavyJetmass_udsc, _h_heavyJetmass_bottom;
    Histo1DPtr _h_totalJetbroad_udsc, _h_totalJetbroad_bottom;
    Histo1DPtr _h_wideJetbroad_udsc, _h_wideJetbroad_bottom;
    Histo1DPtr _h_Cparameter_udsc, _h_Cparameter_bottom;
    Histo1DPtr _h_Dparameter_udsc, _h_Dparameter_bottom;
    Histo1DPtr _h_Ncharged, _h_Ncharged_udsc, _h_Ncharged_bottom;
    Histo1DPtr _h_scaledMomentum, _h_scaledMomentum_udsc, _h_scaledMomentum_bottom;
    // at other enegies
    Histo1DPtr _h_thrust,_h_rho,_h_B_T,_h_B_W,_h_C,_h_D,_h_N,_h_xi;
    // todo add the jets
    // Histo1DPtr _h_y_2_JADE,_h_y_3_JADE,_h_y_4_JADE,_h_y_5_JADE;
    // Histo1DPtr _h_y_2_Durham,_h_y_3_Durham,_h_y_4_Durham,_h_y_5_Durham;
    // Histo1DPtr _h_y_2_Cambridge,_h_y_3_Cambridge,_h_y_4_Cambridge,_h_y_5_Cambridge;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(L3_2004_I652683);

}
