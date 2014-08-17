// -*- C++ -*-
#ifndef RIVET_ParticleFinder_HH
#define RIVET_ParticleFinder_HH

#include "Rivet/Projection.hh"

namespace Rivet {


  /// @brief Base class for projections which return subsets of an event's particles
  class ParticleFinder : public Projection {
  public:

    /// @name Object lifetime management
    //@{

    /// Virtual destructor for inheritance
    virtual ~ParticleFinder() {}

    /// Clone on the heap.
    virtual const Projection* clone() const = 0;

    //@}


    /// @name Particle accessors
    //@{

    /// Access the projected final-state particles.
    virtual size_t size() const { return _theParticles.size(); }

    /// Is this final state empty?
    virtual bool empty() const { return _theParticles.empty(); }
    /// @deprecated Is this final state empty?
    virtual bool isEmpty() const { return _theParticles.empty(); }

    /// Get the final-state particles in no particular order, with no cuts.
    virtual const Particles& particles() const { return _theParticles; }

    /// @brief Get the final-state particles, with optional cuts.
    /// @note Returns a copy rather than a reference, due to cuts
    virtual const Particles particles(double ptmin, double ptmax=MAXDOUBLE,
                                      double rapmin=-MAXDOUBLE, double rapmax=MAXDOUBLE,
                                      RapScheme rapscheme=PSEUDORAPIDITY) const { return _theParticles; }

    /// Get the final-state particles, ordered by supplied sorting function object.
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    /// @todo Update to use Cuts
    template <typename F>
    const Particles particles(F sorter,
                              double ptmin=0.0, double ptmax=MAXDOUBLE,
                              double rapmin=-MAXDOUBLE, double rapmax=MAXDOUBLE,
                              RapScheme rapscheme=PSEUDORAPIDITY) const {
      std::sort(_theParticles.begin(), _theParticles.end(), sorter);
      return _theParticles;
    }

    /// Get the final-state particles, ordered by decreasing \f$ p_T \f$.
    /// @todo Update to use Cuts
    const Particles particlesByPt(double ptmin=0.0, double ptmax=MAXDOUBLE,
                                  double rapmin=-MAXDOUBLE, double rapmax=MAXDOUBLE,
                                  RapScheme rapscheme=PSEUDORAPIDITY) const {
      return particles(cmpMomByPt);
    }

    /// Get the final-state particles, ordered by decreasing \f$ p \f$.
    /// @todo Update to use Cuts
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    const Particles particlesByP(double ptmin=0.0, double ptmax=MAXDOUBLE,
                                 double rapmin=-MAXDOUBLE, double rapmax=MAXDOUBLE,
                                 RapScheme rapscheme=PSEUDORAPIDITY) const {
      return particles(cmpMomByP);
    }

    /// Get the final-state particles, ordered by decreasing \f$ E \f$.
    /// @todo Update to use Cuts
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    const Particles particlesByE(double ptmin=0.0, double ptmax=MAXDOUBLE,
                                 double rapmin=-MAXDOUBLE, double rapmax=MAXDOUBLE,
                                 RapScheme rapscheme=PSEUDORAPIDITY) const {
      return particles(cmpMomByE);
    }

    /// Get the final-state particles, ordered by decreasing \f$ E_T \f$.
    /// @todo Update to use Cuts
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    const Particles particlesByEt(double ptmin=0.0, double ptmax=MAXDOUBLE,
                                  double rapmin=-MAXDOUBLE, double rapmax=MAXDOUBLE,
                                  RapScheme rapscheme=PSEUDORAPIDITY) const {
      return particles(cmpMomByEt);
    }

    /// Get the final-state particles, ordered by increasing \f$ \eta \f$.
    /// @todo Update to use Cuts
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    const Particles particlesByEta(double ptmin=0.0, double ptmax=MAXDOUBLE,
                                   double rapmin=-MAXDOUBLE, double rapmax=MAXDOUBLE,
                                   RapScheme rapscheme=PSEUDORAPIDITY) const {
      return particles(cmpMomByAscPseudorapidity);
    }

    /// Get the final-state particles, ordered by increasing \f$ |\eta| \f$.
    /// @todo Update to use Cuts
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    const Particles particlesByModEta(double ptmin=0.0, double ptmax=MAXDOUBLE,
                                      double rapmin=-MAXDOUBLE, double rapmax=MAXDOUBLE,
                                      RapScheme rapscheme=PSEUDORAPIDITY) const {
      return particles(cmpMomByAscAbsPseudorapidity);
    }

    /// Get the final-state particles, ordered by increasing \f$ y \f$.
    /// @todo Update to use Cuts
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    const Particles particlesByRapidity(double ptmin=0.0, double ptmax=MAXDOUBLE,
                                        double rapmin=-MAXDOUBLE, double rapmax=MAXDOUBLE,
                                        RapScheme rapscheme=PSEUDORAPIDITY) const {
      return particles(cmpMomByAscRapidity);
    }

    /// Get the final-state particles, ordered by increasing \f$ |y| \f$.
    /// @todo Update to use Cuts
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    const Particles particlesByModRapidity(double ptmin=0.0, double ptmax=MAXDOUBLE,
                                           double rapmin=-MAXDOUBLE, double rapmax=MAXDOUBLE,
                                           RapScheme rapscheme=PSEUDORAPIDITY) const {
      return particles(cmpMomByAscAbsRapidity);
    }

    //@}


    /// Minimum-\f$ p_\perp \f$ requirement.
    /// @todo Replace with cuts() accessor
    virtual double ptMin() const { return _ptmin; }


    /// @name JetAlg compatibility
    //@{

    typedef Particle entity_type;
    typedef Particles collection_type;

    /// Template-usable interface common to JetAlg.
    const collection_type& entities() const {
      return particles();
    }

    //@}


  protected:

    /// Apply the projection to the event.
    virtual void project(const Event& e) = 0;

    /// Compare projections.
    virtual int compare(const Projection& p) const;

    /// Decide if a particle is to be accepted or not.
    // virtual bool accept(const Particle& p) const;


  protected:

    /// The ranges allowed for pseudorapidity.
    /// @todo Replace with a single Cuts object
    vector<pair<double,double> > _etaRanges;
    /// The minimum allowed transverse momentum.
    /// @todo Replace with a single Cuts object
    double _ptmin;

    /// The final-state particles.
    mutable Particles _theParticles;

  };


}

#endif
