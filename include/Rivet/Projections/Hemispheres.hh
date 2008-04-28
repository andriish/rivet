// -*- C++ -*-
#ifndef RIVET_Hemispheres_HH
#define RIVET_Hemispheres_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/AxesDefinition.hh"
#include "Rivet/Event.hh"


namespace Rivet {

  /**
     @brief Calculate the hemisphere masses and broadenings.

     Calculate the hemisphere masses and broadenings, with event hemispheres
     defined by the plane normal to the thrust vector, \f$ \vec{n}_\mathrm{T} \f$.

     @todo Allow axes to be defined by sphericity: superclass Thrust and Sphericity as AxisDefinition?

     The "high" hemisphere mass, 
     \f$ M^2_\mathrm{high} / E^2_\mathrm{vis} \f$, is defined as
     \f[
     \frac{M^2_\mathrm{high}}{E^2_\mathrm{vis}} = 
     \frac{1}{E^2_\mathrm{vis}} \max
     \left(
     \left| \sum_{\vec{p}_k \cdot \vec{n}_\mathrm{T} > 0} p_k \right|^2 ,
     \left| \sum_{\vec{p}_k \cdot \vec{n}_\mathrm{T} < 0} p_k \right|^2
     \right)
     \f]
     and the corresponding "low" hemisphere mass, 
     \f$ M^2_\mathrm{low} / E^2_\mathrm{vis} \f$,
     is the sum of momentum vectors in the opposite hemisphere, i.e. 
     \f$ \max \rightarrow \min \f$ in the formula above.

     Finally, we define a hemisphere mass difference:
     \f[
     \frac{M^2_\mathrm{diff} }{ E^2_\mathrm{vis}} = 
     \frac{ M^2_\mathrm{high} - M^2_\mathrm{low} }{ E^2_\mathrm{vis}} .
     \f]

     Similarly to the masses, we also define hemisphere broadenings, using the
     momenta transverse to the thrust axis:
     \f[
     B_\pm =
     \frac{
       \sum{\pm \vec{p}_i \cdot \vec{n}_\mathrm{T} > 0} 
       |\vec{p}_i \times \vec{n}_\mathrm{T} | 
     }{
       2 \sum_i | \vec{p}_i | 
     }
     \f]
     and then a set of the broadening maximum, minimum, sum and difference as follows:
     \f[ B_\mathrm{max}  = \max(B_+, B_-) \f]
     \f[ B_\mathrm{min}  = \min(B_+, B_-) \f]
     \f[ B_\mathrm{sum}  = B_+ + B_- \f]
     \f[ B_\mathrm{diff} = |B_+ - B_-| \f]

     Internally, this projection uses the Thrust projection to determine the
     hemisphere orientation.
  */
  class Hemispheres : public Projection {
  public:

    /// Constructor.
    Hemispheres(const AxesDefinition& ax)
      : _E2vis(-1), _M2high(-1), _M2low(-1),
        _Bmax(-1), _Bmin(-1),
        _highMassEqMaxBroad(true)
    {
      setName("Hemispheres");
      addProjection(ax, "Axes");
    }


  protected:

    /// Perform the projection on the Event.
    void project(const Event& e);

    /// Compare with other projections.
    int compare(const Projection& p) const {
      return mkNamedPCmp(p, "Axes");
    }


  public:

    /// @name Hemisphere masses (scaled by \f$ 1 / E^2_\mathrm{vis} \f$).
    ///@{
    const double getE2vis() const { return _E2vis; }
    const double getM2high() const { return _M2high; }
    const double getM2low() const { return _M2low; }
    const double getM2diff() const { return _M2high -_M2low; }
    const double getScaledM2high() const { 
      if (_M2high == 0.0) return 0.0;
      if (_E2vis != 0.0) return _M2high/_E2vis; 
      else return std::numeric_limits<double>::max(); 
    }
    const double getScaledM2low() const {
      if (_M2low == 0.0) return 0.0;
      if (_E2vis != 0.0) return _M2low/_E2vis;
      else return std::numeric_limits<double>::max(); 
    }
    const double getScaledM2diff() const { 
      if (getM2diff() == 0.0) return 0.0;
      if (_E2vis != 0.0) return getM2diff()/_E2vis; 
      else return std::numeric_limits<double>::max(); 
    }
    ///@}


    /// @name Hemisphere broadenings.
    ///@{
    const double getBmax() const { return _Bmax; }
    const double getBmin() const { return _Bmin; }
    const double getBsum() const { return _Bmax + _Bmin; }
    const double getBdiff() const { return fabs(_Bmax - _Bmin); } // <- fabs(), just in case...
    ///@}


    /// Is the hemisphere with the max mass the same as the one with the max broadening?
    const bool massMatchesBroadening() {
      return _highMassEqMaxBroad;
    }

        
  private:

    /// Visible energy-squared, \f$ E^2_\mathrm{vis} \f$.
    double _E2vis;

    /// Hemisphere mass variables.
    double _M2high, _M2low;

    /// Hemisphere broadening variables.
    double _Bmax, _Bmin;

    /// Is the hemisphere with the max mass the same as the one with the max broadening?
    bool _highMassEqMaxBroad;

  };

}

#endif
