// -*- C++ -*-
#ifndef RIVET_Jet_HH
#define RIVET_Jet_HH

#include "Rivet/Rivet.hh"
#include <numeric>

namespace Rivet {


  /// A minimal class representing a jet of particles.
  class Jet : public ParticleBase {
  public:

    /// Constructor.
    Jet();

    /// Define a Jet::iterator via a typedef.
    typedef vector<FourMomentum>::iterator iterator;

    /// Define a Jet::const_iterator via a typedef.
    typedef vector<FourMomentum>::const_iterator const_iterator;

    /// Get a begin iterator over the particle/track four-momenta in this jet.
    iterator begin() {
      return _particles.begin();
    }

    /// Get an end iterator over the particle/track four-momenta in this jet.
    iterator end() {
      return _particles.end();
    }

    /// Get a const begin iterator over the particle/track four-momenta in this jet.
    const_iterator begin() const {
      return _particles.begin();
    }

    /// Get a const end iterator over the particle/track four-momenta in this jet.
    const_iterator end() const {
      return _particles.end();
    }

    /// Get the track momenta in this jet.
    vector<FourMomentum>& momenta() {
      return _particles;
    }

    /// Get the track momenta in this jet (const version).
    const vector<FourMomentum>& momenta() const {
      return _particles;
    }

    /// Get the Rivet::Particles (full information) in this jet
    vector<Particle>& particles() {
      return _fullParticles;
    }
    
    /// Get the Rivet::Particles (full information) in this jet (const version)
    const vector<Particle>& particles() const {
      return _fullParticles;
    }
    
    /// Number of particles (tracks) in this jet.
    size_t size() const {
      return _particles.size();
    }

    /// Set the particles/tracks collection.
    Jet& setParticles(vector<FourMomentum> particles);

    /// Add a particle/track to this jet.
    Jet& addParticle(FourMomentum particle);

    /// Add a particle/track to this jet.
    Jet& addParticle(const Particle& particle);
    
    /// Check whether this jet contains a particular particle.
    bool containsParticle(const Particle& particle) const;

    /// Check whether this jet contains a certain particle type.
    bool containsParticleId(PdgId pid) const;
 
    /// Reset this jet as empty.
    Jet& clear();

    /// Get the average \f$ \eta \f$ for this jet, with the average weighted
    /// by the \f$ p_T \f$ values of the constituent tracks. (caches)
    double ptWeightedEta() const;

    /// Get the average \f$ \phi \f$ for this jet, with the average weighted
    /// by the \f$ p_T \f$ values of the constituent tracks. (caches)
    double ptWeightedPhi() const;

    /// Get the unweighted average \f$ \eta \f$ for this jet. (caches)
    double eta() const;

    /// Get the unweighted average \f$ \phi \f$ for this jet. (caches)
    double phi() const;

    /// Get equivalent single momentum four-vector. (caches)
    const FourMomentum& momentum() const;
    
    /// Get equivalent single momentum four-vector. (caches)
    FourMomentum& momentum();


  public:
    
    /// Get the sum of the \f$ p_T \f$ values of the constituent tracks. (caches)
    double ptSum() const;

    /// Get the sum of the \f$ E_T \f$ values of the constituent tracks. (caches)
    double EtSum() const;


  private:

    /// Clear the internal cached values.
    void _resetCaches() const;

    // Calculate cached equivalent momentum vector
    void _calcMomVector() const;

    /// Internal caching method to calculate the average \f$ \eta \f$ and \f$
    /// \phi \f$ for this jet, weighted by the \f$ p_T \f$ values of the
    /// constituent tracks.
    void _calcPtAvgs() const;

    /// Internal caching method to calculate the unweighted average \f$ \eta
    /// \f$ and \f$ \phi \f$ for this jet.
    void _calcAvgs() const;


  private:

    /// The particle tracks.
    std::vector<FourMomentum> _particles;

    /// Full particle information including tracks, ID etc
    ParticleVector _fullParticles;
    
    /// Cached values of \f$ \bar{\phi} \f$ and \f$ \bar{\eta} \f$.
    mutable double _phi, _eta;
    mutable bool _okPhi, _okEta;

    /// Cached values of the \f$ p_T \f$-weighted \f$ \bar{\phi} \f$ and \f$ \bar{\eta} \f$.
    mutable double _ptWeightedPhi, _ptWeightedEta;
    mutable bool _okPtWeightedPhi, _okPtWeightedEta;

    /// Cached value of the \f$ p_T \f$ sum.
    mutable double _totalPt;
    mutable bool _okTotalPt;

    /// Cached value of the \f$ E_T \f$ sum.
    mutable double _totalEt;
    mutable bool _okTotalEt;
    
    mutable FourMomentum _momentum;
    mutable bool _okMomentum;
    
  };


  /// Typedef for a collection of Jet objects.
  typedef std::vector<Jet> Jets;

}

#endif
