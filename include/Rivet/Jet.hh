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

    /// Check whether this jet contains a charm-flavoured hadron.
    bool containsCharm() const;

    /// Check whether this jet contains a bottom-flavoured hadron.
    bool containsBottom() const;
 
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
    /// @todo Review if these caches are needed/consistent: just the vector, maybe?
    void _calcPtAvgs() const;

    /// Internal caching method to calculate the unweighted average \f$ \eta
    /// \f$ and \f$ \phi \f$ for this jet.
    /// @todo Review if these caches are needed/consistent: just the vector, maybe?
    void _calcAvgs() const;


  private:

    /// The particle tracks.
    std::vector<FourMomentum> _particles;

    /// Full particle information including tracks, ID etc
    ParticleVector _fullParticles;
    
    /// Cached values of \f$ \bar{\phi} \f$ and \f$ \bar{\eta} \f$.
    /// @todo Review if these caches are needed/consistent: just the vector, maybe?
    mutable double _phi, _eta;
    /// @todo Review if these caches are needed/consistent: just the vector, maybe?
    mutable bool _okPhi, _okEta;

    /// Cached values of the \f$ p_T \f$-weighted \f$ \bar{\phi} \f$ and \f$ \bar{\eta} \f$.
    /// @todo Review if these caches are needed/consistent: just the vector, maybe?
    mutable double _ptWeightedPhi, _ptWeightedEta;
    /// @todo Review if these caches are needed/consistent: just the vector, maybe?
    mutable bool _okPtWeightedPhi, _okPtWeightedEta;

    /// Cached value of the \f$ p_T \f$ sum.
    /// @todo Review if these caches are needed/consistent: just the vector, maybe?
    mutable double _totalPt;
    /// @todo Review if these caches are needed/consistent: just the vector, maybe?
    mutable bool _okTotalPt;

    /// Cached value of the \f$ E_T \f$ sum.
    /// @todo Review if these caches are needed/consistent: just the vector, maybe?
    mutable double _totalEt;
    /// @todo Review if these caches are needed/consistent: just the vector, maybe?
    mutable bool _okTotalEt;
    
    mutable FourMomentum _momentum;
    mutable bool _okMomentum;
    
  };


  /// Typedef for a collection of Jet objects.
  typedef std::vector<Jet> Jets;


  /// @name Jet comparison functions for STL sorting
  //@{

  // @brief Compare jets by \f$ p_\perp \f$ (descending - usual sorting for HEP)
  /// Use this so that highest \f$ p_\perp \f$ is at the front of the list
  inline bool cmpJetsByPt(const Jet& a, const Jet& b) {
    return a.ptSum() > b.ptSum();
  }
  // @brief Compare jets by \f$ p_\perp \f$ (ascending)
  /// Use this so that lowest \f$ p_\perp \f$ is at the front of the list
  inline bool cmpJetsByAscPt(const Jet& a, const Jet& b) {
    return a.ptSum() < b.ptSum();
  }

  // @brief Compare jets by \f$ E_\perp \f$ (descending - usual sorting for HEP)
  /// Use this so that highest \f$ E_\perp \f$ is at the front of the list
  inline bool cmpJetsByEt(const Jet& a, const Jet& b) {
    return a.EtSum() > b.EtSum();
  }
  // @brief Compare jets by \f$ E_\perp \f$ (ascending)
  /// Use this so that lowest \f$ E_\perp \f$ is at the front of the list
  inline bool cmpJetsByEtDesc(const Jet& a, const Jet& b) {
    return a.EtSum() < b.EtSum();
  }

  // @brief Compare jets by \f$ E \f$ (descending - usual sorting for HEP)
  /// Use this so that highest \f$ E \f$ is at the front of the list
  inline bool cmpJetsByE(const Jet& a, const Jet& b) {
    return a.momentum().E() > b.momentum().E();
  }
  // @brief Compare jets by \f$ E \f$ (ascending)
  /// Use this so that lowest \f$ E \f$ is at the front of the list
  inline bool cmpJetsByAscE(const Jet& a, const Jet& b) {
    return a.momentum().E() < b.momentum().E();
  }

  // @brief Compare jets by \f$ \eta \f$ (descending)
  /// Use this so that highest \f$ \eta \f$ is at the front of the list
  inline bool cmpJetsByDescPseudorapidity(const Jet& a, const Jet& b) {
    return a.momentum().pseudorapidity() > b.momentum().pseudorapidity();
  }
  // @brief Compare jets by \f$ \eta \f$ (ascending)
  /// Use this so that lowest \f$ \eta \f$ is at the front of the list
  inline bool cmpJetsByAscPseudorapidity(const Jet& a, const Jet& b) {
    return a.momentum().pseudorapidity() < b.momentum().pseudorapidity();
  }

  // @brief Compare jets by \f$ |\eta| \f$ (descending)
  /// Use this so that highest \f$ |\eta| \f$ is at the front of the list
  inline bool cmpJetsByDescAbsPseudorapidity(const Jet& a, const Jet& b) {
    return fabs(a.momentum().pseudorapidity()) > fabs(b.momentum().pseudorapidity());
  }
  // @brief Compare jets by \f$ |\eta| \f$ (ascending)
  /// Use this so that lowest \f$ |\eta| \f$ is at the front of the list
  inline bool cmpJetsByAscAbsPseudorapidity(const Jet& a, const Jet& b) {
    return fabs(a.momentum().pseudorapidity()) < fabs(b.momentum().pseudorapidity());
  }

  // @brief Compare jets by \f$ y \f$ (descending)
  /// Use this so that highest \f$ y \f$ is at the front of the list
  inline bool cmpJetsByDescRapidity(const Jet& a, const Jet& b) {
    return a.momentum().rapidity() > b.momentum().rapidity();
  }
  // @brief Compare jets by \f$ y \f$ (ascending)
  /// Use this so that lowest \f$ y \f$ is at the front of the list
  inline bool cmpJetsByAscRapidity(const Jet& a, const Jet& b) {
    return a.momentum().rapidity() < b.momentum().rapidity();
  }

  // @brief Compare jets by \f$ |y| \f$ (descending)
  /// Use this so that highest \f$ |y| \f$ is at the front of the list
  inline bool cmpJetsByDescAbsRapidity(const Jet& a, const Jet& b) {
    return fabs(a.momentum().rapidity()) > fabs(b.momentum().rapidity());
  }
  // @brief Compare jets by \f$ |y| \f$ (ascending)
  /// Use this so that lowest \f$ |y| \f$ is at the front of the list
  inline bool cmpJetsByAscAbsRapidity(const Jet& a, const Jet& b) {
    return fabs(a.momentum().rapidity()) < fabs(b.momentum().rapidity());
  }

  //@}

}

#endif
