// -*- C++ -*-
#ifndef RIVET_GammaGammaFinalState_HH
#define RIVET_GammaGammaFinalState_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/GammaGammaKinematics.hh"

namespace Rivet {


  /// @brief Final state particles boosted to the hadronic center of mass system.
  ///
  /// NB. The GammaGamma scattered lepton is not included in the final state particles.
  class GammaGammaFinalState: public FinalState {
  public:

    /// Type of GammaGamma boost to apply
    enum BoostType { HCM, BREIT, LAB };


    /// @name Constructors
    //@{

    /// Constructor with explicit FinalState
    /// @note The GammaGammaKinematics has no parameters, hence explicitly passing it as an arg shouldn't be necessary.
    GammaGammaFinalState(const FinalState& fs, BoostType boosttype, const GammaGammaKinematics& kinematicsp=GammaGammaKinematics())
      : _boosttype(boosttype)
    {
      setName("GammaGammaFinalState");
      declare(fs, "FS");
      declare(kinematicsp, "Kinematics");
    }

    /// Constructor with optional FinalState
    /// @note The GammaGammaKinematics has no parameters, hence explicitly passing it as an arg shouldn't be necessary.
    GammaGammaFinalState(BoostType boosttype, const FinalState& fs=FinalState(), const GammaGammaKinematics& kinematicsp=GammaGammaKinematics())
      : GammaGammaFinalState(fs, boosttype, kinematicsp)
    {    }

    /// Constructor with explicit cuts to define final-state particles
    /// @note The GammaGammaKinematics has no parameters, hence explicitly passing it as an arg shouldn't be necessary.
    GammaGammaFinalState(const Cut& c, BoostType boosttype, const GammaGammaKinematics& kinematicsp=GammaGammaKinematics())
      : GammaGammaFinalState(FinalState(c), boosttype, kinematicsp)
    {    }

    /// Constructor with explicit cuts to define final-state particles
    /// @note The GammaGammaKinematics has no parameters, hence explicitly passing it as an arg shouldn't be necessary.
    GammaGammaFinalState(BoostType boosttype, const Cut& c, const GammaGammaKinematics& kinematicsp=GammaGammaKinematics())
      : GammaGammaFinalState(FinalState(c), boosttype, kinematicsp)
    {    }

    // /// @brief Constructor with default FinalState
    // /// @note The GammaGammaKinematics has no parameters, hence explicitly passing it as an arg shouldn't be necessary.
    // GammaGammaFinalState(BoostType boosttype, const GammaGammaKinematics& kinematicsp=GammaGammaKinematics())
    //   : GammaGammaFinalState(FinalState(), boosttype, kinematicsp)
    // {    }

    /// Backward compatible constructor with default FinalState
    /// @deprecated Prefer a version that doesn't need a GammaGammaKinematics argument
    GammaGammaFinalState(const GammaGammaKinematics& kinematicsp, BoostType boosttype)
      : GammaGammaFinalState(FinalState(), boosttype, kinematicsp)
    {    }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(GammaGammaFinalState);

    //@}


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const {
      const GammaGammaFinalState& other = dynamic_cast<const GammaGammaFinalState&>(p);
      return mkNamedPCmp(p, "Kinematics") || mkNamedPCmp(p, "FS") || cmp(_boosttype, other._boosttype);
    }


  private:

    BoostType _boosttype;

  };


}

#endif
