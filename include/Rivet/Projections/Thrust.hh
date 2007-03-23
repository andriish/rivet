// -*- C++ -*-
#ifndef RIVET_Thrust_HH
#define RIVET_Thrust_HH
// Declaration of the Thrust class.

#include "Rivet/Projections/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
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
    inline Thrust(FinalState& fsp) 
      : _calculatedThrust(false), _fsproj(&fsp)
    { }

    /// Copy constructor.
    /// @todo More? or is the default enough?
    inline Thrust(const Thrust& x)
      : Projection(x)
    { }

    /// Destructor.
    virtual ~Thrust() { }
    //@}

  public:
    /// Return the name of the projection
    inline string name() const {
      return "Thrust";
    }

  protected:

    /// Perform the projection on the Event: only to be called by 
    /// Event::applyProjection(Projection &).
    void project(const Event& e);

    /// This function defines a unique ordering between different 
    /// Projection objects of the same class. Should only be called from 
    /// operator<(const Projection &).
    int compare(const Projection& p) const;

  public:

    /// The thrust scalar (maximum thrust).
    inline const double thrust() const { return _thrusts[0]; }

    /// The thrust major scalar (thrust along thrust major axis).
    inline const double thrustMajor() const { return _thrusts[1]; }

    /// The thrust minor scalar (thrust along thrust minor axis).
    inline const double thrustMinor() const { return _thrusts[2]; }

    /// The thrust axis.
    inline const Vector3& thrustAxis() const { return _thrustAxes[0]; }

    /// The thrust major axis (axis of max thrust perpendicular to thrust axis).
    inline const Vector3& thrustMajorAxis() const { return _thrustAxes[1]; }

    /// The thrust minor axis (axis perpendicular to thrust and thrust major).
    inline const Vector3& thrustMinorAxis() const { return _thrustAxes[2]; }


  private:

    /// The thrust scalars.
    vector<double> _thrusts;

    /// The thrust axes.
    vector<Vector3> _thrustAxes;

    /// Caching flag to avoid costly recalculations.
    bool _calculatedThrust;

    /// The FinalState projection used by this projection
    FinalState* _fsproj;

  private:

    /// Calculate the thrust axes and scalars, using special cases where possible.
    void calcThrust(const FinalState& fs);

    /// Explicitly calculate the thrust minor value and axis.
    void calcM(const vector<Vector3>& p, double& m, Vector3& maxis) const;

    /// Explicitly calculate the thrust value and axis.
    void calcT(const vector<Vector3>& p, double& t, Vector3& taxis) const;

    /// The assignment operator is private and must never be called.
    Thrust& operator=(const Thrust &);

  };
  
}


#endif /* RIVET_Thrust_HH */
