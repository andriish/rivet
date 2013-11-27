#ifndef RIVET_ParticleBase_HH
#define RIVET_ParticleBase_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Math/Vectors.hh"

namespace Rivet {


  /// @brief Base class for particle-like things like Particle and Jet
  class ParticleBase {
  public:

    /// Default constructor
    ParticleBase() { }

    /// Virtual destructor
    virtual ~ParticleBase() { }


    /// @name Effective momentum accessors
    //@{

    /// Get equivalent single momentum four-vector.
    // virtual FourMomentum& momentum() = 0;

    /// Get equivalent single momentum four-vector (const).
    virtual const FourMomentum& momentum() const = 0;

    //@}


    /// @name Convenience access to the effective 4-vector properties
    //@{

    /// Get the energy directly.
    double energy() const { return momentum().E(); }
    /// Get the energy directly.
    double E() const { return momentum().E(); }

    /// Get the \f$ p_T \f$ directly.
    double pT() const { return momentum().pT(); }
    /// Get the \f$ p_T \f$ directly (alias).
    double perp() const { return momentum().pT(); }

    /// Get the \f$ E_T \f$ directly.
    double Et() const { return momentum().Et(); }

    /// Get the mass directly.
    double mass() const { return momentum().mass(); }
    /// Get the mass**2 directly.
    double mass2() const { return momentum().mass2(); }

    /// Get the \f$ \eta \f$ directly.
    double pseudorapidity() const { return momentum().eta(); }
    /// Get the \f$ \eta \f$ directly (alias).
    double eta() const { return momentum().eta(); }
    /// Get the \f$ |\eta| \f$ directly.
    double abspseudorapidity() const { return fabs(momentum().eta()); }
    /// Get the \f$ |\eta| \f$ directly (alias).
    double abseta() const { return fabs(momentum().eta()); }

    /// Get the \f$ y \f$ directly.
    double rapidity() const { return momentum().rapidity(); }
    /// Get the \f$ y \f$ directly (alias).
    double rap() const { return momentum().rapidity(); }
    /// Get the \f$ |y| \f$ directly.
    double absrapidity() const { return fabs(momentum().rapidity()); }
    /// Get the \f$ |y| \f$ directly (alias).
    double absrap() const { return fabs(momentum().rapidity()); }

    /// Get the \f$ \phi \f$ directly.
    double phi() const { return momentum().phi(); }

    //@}



    /// Struct for sorting by increasing transverse momentum in STL set, sort, etc.
    struct byPTAscending {
      bool operator()(const ParticleBase& left, const ParticleBase& right) const {
        double pt2left = left.momentum().pT2();
        double pt2right = right.momentum().pT2();
        return pt2left < pt2right;
      }

      bool operator()(const ParticleBase* left, const ParticleBase* right) const {
        return (*this)(*left, *right);
      }
    };


    /// Struct for sorting by decreasing transverse momentum in STL set, sort etc.
    struct byPTDescending {
      bool operator()(const ParticleBase& left, const ParticleBase& right) const {
        return byPTAscending()(right, left);
      }

      bool operator()(const ParticleBase* left, const ParticleBase* right) const {
        return (*this)(*left, *right);
      }
    };


    /// Struct for sorting by increasing transverse energy
    struct byETAscending {
      bool operator()(const ParticleBase& left, const ParticleBase& right) const {
        double pt2left = left.momentum().Et2();
        double pt2right = right.momentum().Et2();
        return pt2left < pt2right;
      }

      bool operator()(const ParticleBase* left, const ParticleBase* right) const {
        return (*this)(*left, *right);
      }
    };


    /// Struct for sorting by decreasing transverse energy
    struct byETDescending {
      bool operator()(const ParticleBase& left, const ParticleBase& right) const {
        return byETAscending()(right, left);
      }

      bool operator()(const ParticleBase* left, const ParticleBase* right) const {
        return (*this)(*left, *right);
      }
    };


    /// Struct for sorting by increasing energy
    struct byEAscending {
      bool operator()(const ParticleBase& left, const ParticleBase& right) const {
        double pt2left = left.momentum().E();
        double pt2right = right.momentum().E();
        return pt2left < pt2right;
      }

      bool operator()(const ParticleBase* left, const ParticleBase* right) const {
        return (*this)(*left, *right);
      }
    };


    /// Struct for sorting by decreasing energy
    struct byEDescending {
      bool operator()(const ParticleBase& left, const ParticleBase& right) const {
        return byEAscending()(right, left);
      }

      bool operator()(const ParticleBase* left, const ParticleBase* right) const {
        return (*this)(*left, *right);
      }
    };


  };


}

#endif
