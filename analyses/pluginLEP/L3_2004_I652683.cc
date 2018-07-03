// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"

#define I_KNOW_THE_INITIAL_QUARKS_PROJECTION_IS_DODGY_BUT_NEED_TO_USE_IT
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  /// Jet rates and event shapes at LEP I+II
  class L3_2004_I652683 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(L3_2004_I652683);
    // L3_2004_I652683() : Analysis("L3_2004_I652683")
    // {    }


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
      book(_h_Thrust_udsc          , 47, 1, 1);
      book(_h_Thrust_bottom        , 47, 1, 2);
      book(_h_heavyJetmass_udsc    , 48, 1, 1);
      book(_h_heavyJetmass_bottom  , 48, 1, 2);
      book(_h_totalJetbroad_udsc   , 49, 1, 1);
      book(_h_totalJetbroad_bottom , 49, 1, 2);
      book(_h_wideJetbroad_udsc    , 50, 1, 1);
      book(_h_wideJetbroad_bottom  , 50, 1, 2);
      book(_h_Cparameter_udsc      , 51, 1, 1);
      book(_h_Cparameter_bottom    , 51, 1, 2);
      book(_h_Dparameter_udsc      , 52, 1, 1);
      book(_h_Dparameter_bottom    , 52, 1, 2);
      book(_h_Ncharged             , 59, 1, 1);
      book(_h_Ncharged_udsc        , 59, 1, 2);
      book(_h_Ncharged_bottom      , 59, 1, 3);
      book(_h_scaledMomentum       , 65, 1, 1);
      book(_h_scaledMomentum_udsc  , 65, 1, 2);
      book(_h_scaledMomentum_bottom, 65, 1, 3);

      book(_sumW_udsc, "sumW_udsc");
      book(_sumW_b, "sumW_b");
      book(_sumW_ch, "sumW_ch");
      book(_sumW_ch_udsc, "sumW_ch_udsc");
      book(_sumW_ch_b, "sumW_ch_b");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get beam average momentum
      const ParticlePair& beams = apply<Beam>(event, "beams").beams();
      const double beamMomentum = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;

      // InitialQuarks projection to have udsc events separated from b events
      /// @todo Yuck!!! Eliminate when possible...
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
      const int iflav = (flavour == PID::DQUARK || flavour == PID::UQUARK || flavour == PID::SQUARK || flavour == PID::CQUARK) ? 1 : (flavour == PID::BQUARK) ? 5 : 0;

      // Update weight sums
      if (iflav == 1) {
        _sumW_udsc->fill();
      } else if (iflav == 5) {
        _sumW_b->fill();
      }
      _sumW_ch->fill();

      // Charged multiplicity
      const FinalState& cfs = applyProjection<FinalState>(event, "CFS");
      _h_Ncharged->fill(cfs.size());
      if (iflav == 1) {
        _sumW_ch_udsc->fill();
        _h_Ncharged_udsc->fill(cfs.size());
      } else if (iflav == 5) {
        _sumW_ch_b->fill();
        _h_Ncharged_bottom->fill(cfs.size());
      }

      // Scaled momentum
      const Particles& chparticles = cfs.particlesByPt();
      for (const Particle& p : chparticles) {
        const Vector3 momentum3 = p.p3();
        const double mom = momentum3.mod();
        const double scaledMom = mom/beamMomentum;
        const double logScaledMom = std::log(scaledMom);
        _h_scaledMomentum->fill(-logScaledMom);
        if (iflav == 1) {
          _h_scaledMomentum_udsc->fill(-logScaledMom);
        } else if (iflav == 5) {
          _h_scaledMomentum_bottom->fill(-logScaledMom);
        }
      }

      // Thrust
      const Thrust& thrust = applyProjection<Thrust>(event, "thrust");
      if (iflav == 1) {
        _h_Thrust_udsc->fill(thrust.thrust());
      } else if (iflav == 5) {
        _h_Thrust_bottom->fill(thrust.thrust());
      }

      // C and D Parisi parameters
      const ParisiTensor& parisi = applyProjection<ParisiTensor>(event, "Parisi");
      if (iflav == 1) {
        _h_Cparameter_udsc->fill(parisi.C());
        _h_Dparameter_udsc->fill(parisi.D());
      } else if (iflav == 5) {
        _h_Cparameter_bottom->fill(parisi.C());
        _h_Dparameter_bottom->fill(parisi.D());
      }

      // The hemisphere variables
      const Hemispheres& hemisphere = applyProjection<Hemispheres>(event, "Hemispheres");
      if (iflav == 1) {
        _h_heavyJetmass_udsc->fill(hemisphere.scaledM2high());
        _h_totalJetbroad_udsc->fill(hemisphere.Bsum());
        _h_wideJetbroad_udsc->fill(hemisphere.Bmax());
      } else if (iflav == 5) {
        _h_heavyJetmass_bottom->fill(hemisphere.scaledM2high());
        _h_totalJetbroad_bottom->fill(hemisphere.Bsum());
        _h_wideJetbroad_bottom->fill(hemisphere.Bmax());
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale({_h_Thrust_udsc, _h_heavyJetmass_udsc, _h_totalJetbroad_udsc, _h_wideJetbroad_udsc, _h_Cparameter_udsc, _h_Dparameter_udsc}, 1/ *_sumW_udsc);
      scale({_h_Thrust_bottom, _h_heavyJetmass_bottom, _h_totalJetbroad_bottom, _h_wideJetbroad_bottom, _h_Cparameter_bottom, _h_Dparameter_bottom}, 1./ *_sumW_b);
      scale(_h_Ncharged, 2/ *_sumW_ch);
      scale(_h_Ncharged_udsc, 2/ *_sumW_ch_udsc);
      scale(_h_Ncharged_bottom, 2/ *_sumW_ch_b);
      scale(_h_scaledMomentum, 1/ *_sumW_ch);
      scale(_h_scaledMomentum_udsc, 1/ *_sumW_ch_udsc);
      scale(_h_scaledMomentum_bottom, 1/ *_sumW_ch_b);
    }


    /// Weight counters
    CounterPtr _sumW_udsc, _sumW_b, _sumW_ch, _sumW_ch_udsc, _sumW_ch_b;

    /// @name Histograms
    //@{
    Histo1DPtr _h_Thrust_udsc, _h_Thrust_bottom;
    Histo1DPtr _h_heavyJetmass_udsc, _h_heavyJetmass_bottom;
    Histo1DPtr _h_totalJetbroad_udsc, _h_totalJetbroad_bottom;
    Histo1DPtr _h_wideJetbroad_udsc, _h_wideJetbroad_bottom;
    Histo1DPtr _h_Cparameter_udsc, _h_Cparameter_bottom;
    Histo1DPtr _h_Dparameter_udsc, _h_Dparameter_bottom;
    Histo1DPtr _h_Ncharged, _h_Ncharged_udsc, _h_Ncharged_bottom;
    Histo1DPtr _h_scaledMomentum, _h_scaledMomentum_udsc, _h_scaledMomentum_bottom;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(L3_2004_I652683);

}
