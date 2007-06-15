// -*- C++ -*-
#ifndef RIVET_Thrust_HH
#define RIVET_Thrust_HH

#include "Rivet/Projections/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Event.hh"
#include "Rivet/RivetCLHEP.hh"


namespace Rivet {

  /**
    @brief Obtain the e+ e- thrust event shape, consisting of the thrust basis and the 
    thrust scalar values in each direction (the thrust, thrust major and thrust
    minor).

    @author Andy Buckley
   
    The scalar (maximum) thrust is defined as
    \f[
    T = \mathrm{max}_{\vec{n}} \frac{\sum_i \left|\vec{p}_i \cdot \vec{n} \right|}{\sum_i |\vec{p}_i|}
    \f],
    with the direction of the unit vector \f$ \vec{n} \f$ which maximises \f$ T \f$ 
    being identified as the thrust axis. The unit vector which maximises the thrust
    scalar in the plane perpendicular to \f$ \vec{n} \f$ is the "thrust major"
    direction, and the vector perpendicular to both the thrust and thrust major directions
    is the thrust minor. Both the major and minor directions have associated thrust 
    scalars.

    Thrust calculations have particularly simple forms for less than 4 particles, and
    in those cases this projection is computationally minimal. For 4 or more particles,
    a more general calculation must be carried out, based on the Brandt/Dahmen method 
    from Z. Phys. C1 (1978). While a polynomial improvement on the exponential scaling
    of the naive method, this algorithm scales asymptotically as 
    \f$ \mathcal{O}\left( n^3 \right) \f$. Be aware that the thrust may easily be the 
    most computationally demanding projection in Rivet for large events!

    The Rivet implementation of thrust is based heavily on Stefan Gieseke's Herwig++
    re-coding of the 'tasso' code from HERWIG.
   */
  class Thrust : public Projection {

  public:

    /// Constructor. The FinalState projection must live throughout the run.
    inline Thrust(FinalState& fsp)
      : _calculatedThrust(false), _fsproj(&fsp)
    { 
      addProjection(fsp);
    }


  public:
    /// Return the name of the projection
    inline string getName() const {
      return "Thrust";
    }

  protected:

    /// Perform the projection on the Event
    void project(const Event& e);

    /// Compare projections
    int compare(const Projection& p) const;

  public:

    ///@{ Thrust scalar accessors
    /// The thrust scalar (maximum thrust).
    inline const double thrust() const { return _thrusts[0]; }
    /// The thrust major scalar (thrust along thrust major axis).
    inline const double thrustMajor() const { return _thrusts[1]; }
    /// The thrust minor scalar (thrust along thrust minor axis).
    inline const double thrustMinor() const { return _thrusts[2]; }
    ///@}

    ///@{ Thrust axis accessors
    /// The thrust axis.
    inline const Vector3& thrustAxis() const { return _thrustAxes[0]; }
    /// The thrust major axis (axis of max thrust perpendicular to thrust axis).
    inline const Vector3& thrustMajorAxis() const { return _thrustAxes[1]; }
    /// The thrust minor axis (axis perpendicular to thrust and thrust major).
    inline const Vector3& thrustMinorAxis() const { return _thrustAxes[2]; }
    ///@}


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


#endif
