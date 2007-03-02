// -*- C++ -*-
#ifndef RIVET_Thrust_HH
#define RIVET_Thrust_HH
// Declaration of the Thrust class.

#include "Rivet/Projections/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Tools/Event/Particle.hh"
#include "Rivet/Tools/Event/Event.hh"


namespace Rivet {

  /**
   * Project out the event shape.
   */
  class Thrust : public Projection {

  public:

    /// @name Standard constructors and destructors.
    //@{
    /// Default constructor. Must specify a FinalState projection which is
    //assumed to live throughout the run.
    inline Thrust(FinalState & fsp);

    /// Copy constructor.
    inline Thrust(const Thrust &);

    /// Destructor.
    virtual ~Thrust();
    //@}

  protected:

    /// Perform the projection on the Event: only to be called by 
    /// Event::applyProjection(Projection &).
    void project(const Event & e);

    /// This function defines a unique ordering between different 
    /// Projection objects of the same class. Should only be called from 
    /// operator<(const Projection &).
    int compare(const Projection & p) const;

  public:

    /// The thrust scalar (maximum thrust)
    inline const double thrust() const;

    /// The thrust major scalar (thrust along thrust major axis)
    inline const double thrustMajor() const;

    /// The thrust minor scalar (thrust along thrust minor axis)
    inline const double thrustMinor() const;

    /// The thrust axis
    inline const Vector3& thrustAxis() const;

    /// The thrust major axis (axis of max thrust perpendicular to thrust axis)
    inline const Vector3& thrustMajorAxis() const;

    /// The thrust minor axis (axis perpendicular to thrust and thrust major)
    inline const Vector3& thrustMinorAxis() const;


  private:

    /// The thrust scalars and axes
    vector<double> thrusts_;
    vector<Vector3> thrustAxes_;

    /// The FinalState projection used by this projection
    FinalState * fsproj;

  private:

    /// The assignment operator is private and must never be called
    Thrust& operator=(const Thrust &);

  };
  
}

#include "Rivet/Projections/Thrust.icc"

#endif /* RIVET_Thrust_HH */
