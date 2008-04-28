// -*- C++ -*-
#ifndef RIVET_CentralEtHCM_HH
#define RIVET_CentralEtHCM_HH

#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projections/FinalStateHCM.hh"


namespace Rivet {


  /// Sum up Et of all particles in the hadronic final state in the
  /// central rapidity bin of the HCM system.
  class CentralEtHCM : public Projection {

  public:

    /// The default constructor. Must specify a FinalStateHCM projection
    /// object which is guaranteed to live throughout the run.
    CentralEtHCM(const FinalStateHCM& fs)
    {
      setName("CentralEtHCM");
      addProjection(fs, "FS"); 
      addCut("eta", MORE_EQ, -0.5);
      addCut("eta", LESS_EQ,  0.5);
    }

  protected:

    /// Apply the projection on to the Event.
    void project(const Event& e);

    /// Compare with other projections
    int compare(const Projection& p) const;

  public:

    /// The sum of the Et in the central rapidity bin.
    double sumEt() const { return _sumet; }

  private:

    /// The sum of the Et in the central rapidity bin.
    double _sumet;

  };

}


#endif
