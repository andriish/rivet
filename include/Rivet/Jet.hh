// -*- C++ -*-
#ifndef RIVET_Jet_HH
#define RIVET_Jet_HH

#include "Rivet/Rivet.hh"
#include <numeric>

namespace Rivet {


  /// A minimal class representing a jet of particles.
  class Jet {
  public:
    /// Constructor.
    Jet() { 
      clear();
    }

    /// Define a Jet::iterator via a typedef.
    typedef vector<FourMomentum>::iterator iterator;

    /// Define a Jet::const_iterator via a typedef.
    typedef vector<FourMomentum>::const_iterator const_iterator;

    /// Get a begin iterator over the particles/tracks in this jet.
    iterator begin() {
      return _particles.begin();
    }

    /// Get an end iterator over the particles/tracks in this jet.
    iterator end() {
      return _particles.end();
    }

    /// Get a const begin iterator over the particles/tracks in this jet.
    const_iterator begin() const {
      return _particles.begin();
    }

    /// Get a const end iterator over the particles/tracks in this jet.
    const_iterator end() const {
      return _particles.end();
    }

    /// Get the particles (tracks) in this jet.
    vector<FourMomentum>& getParticles() {
      return _particles;
    }

    /// Get the particles (tracks) in this jet (const version).
    const vector<FourMomentum>& getParticles() const {
      return _particles;
    }

    /// Number of particles (tracks) in this jet.
    size_t getNumParticles() const {
      return _particles.size();
    }

    /// Set the particles/tracks collection.
    Jet setParticles(vector<FourMomentum> particles) {
      _particles = particles;
      _resetCaches();
      return *this;
    }

    /// Add a particle/track to this jet.
    Jet addParticle(FourMomentum particle) {
      _particles.push_back(particle);
      _resetCaches();
      return *this;
    }

    /// Reset this jet as empty.
    Jet clear() {
      _particles.clear();
      _resetCaches();
      return *this;
    }

    /// Get the average \f$ \eta \f$ for this jet, with the average weighted
    /// by the \f$ p_T \f$ values of the constituent tracks. (caches)
    double getPtWeightedEta() const {
      _calcPtAvgs();
      assert(_okPtWeightedEta);
      return _ptWeightedEta;
    }

    /// Get the average \f$ \phi \f$ for this jet, with the average weighted
    /// by the \f$ p_T \f$ values of the constituent tracks. (caches)
    double getPtWeightedPhi() const {
      _calcPtAvgs();
      assert(_okPtWeightedPhi);
      return _ptWeightedPhi;
    }

    /// Get the unweighted average \f$ \eta \f$ for this jet. (caches)
    double getEta() const {
      _calcAvgs();
      assert(_okEta);
      return _eta;
    }


    /// Get the unweighted average \f$ \phi \f$ for this jet. (caches)
    double getPhi() const {
      _calcAvgs();
      assert(_okPhi);
      return _phi;
    }


    /// Get a FourMomentum vector that characterises this jet (the sum of the
    /// constitutent momenta).
    FourMomentum vector() const {
      return accumulate(begin(), end(), FourMomentum());
    }


    /// Get the sum of the \f$ p_T \f$ values of the constituent tracks. (caches)
    double getPtSum() const {
      if (!_okTotalPt) {
        double ptsum(0.0);
        for (const_iterator p = this->begin(); p != this->end(); ++p) {
          ptsum += p->pT();
        }
        _totalPt = ptsum;
        _okTotalPt = true;
      }
      return _totalPt;
    }


    /// Get the sum of the \f$ p_T \f$ values of the constituent tracks. (caches)
    double getEtSum() const {
      /// @todo Correct this for the real Et definition.
      return getPtSum();
    }


  private:

    /// Clear the internal cached values.
    void _resetCaches() const {
      if (_okPhi || _okEta || _okPtWeightedPhi || _okPtWeightedEta || _okTotalPt) {
        _okPhi = false;
        _okEta = false;
        _okPtWeightedPhi = false;
        _okPtWeightedEta = false;
        _okTotalPt = false;
      }
    }

    /// Internal caching method to calculate the average \f$ \eta \f$ and \f$
    /// \phi \f$ for this jet, weighted by the \f$ p_T \f$ values of the
    /// constituent tracks.
    void _calcPtAvgs() const {
      if (!_okPtWeightedEta || !_okPtWeightedPhi) {
        double ptwetasum(0.0), ptwphisum(0.0), ptsum(0.0);
	double phibegin = (0.0);
        for (const_iterator p = this->begin(); p != this->end(); ++p) {

          double pt = p->pT();
          ptsum += pt;
          ptwetasum += pt * p->pseudorapidity();

	  if (p == this->begin())
	    phibegin = p->azimuthalAngle();
	  else {
	    double dphi = p->azimuthalAngle() - phibegin;
	    ptwphisum += pt*
	      ( dphi > PI ? dphi-TWOPI :
		dphi < -PI ? dphi+TWOPI : dphi);   
	  }
        }
        _totalPt = ptsum;
        _okTotalPt = true;
        _ptWeightedEta = ptwetasum / getPtSum();
        _okPtWeightedEta = true;
	_ptWeightedPhi = phibegin + ptwphisum / getPtSum();
        _okPtWeightedPhi = true;
      }
    }

    /// Internal caching method to calculate the unweighted average \f$ \eta
    /// \f$ and \f$ \phi \f$ for this jet.
    void _calcAvgs() const {
      if (!_okEta || !_okPhi) {
        double etasum(0.0), phisum(0.0);
        for (const_iterator p = this->begin(); p != this->end(); ++p) {
          etasum += p->pseudorapidity();
          phisum += p->azimuthalAngle();
        }
        const double dnum( getNumParticles() );
        _eta = etasum / dnum;
        _okEta = true;
        _phi = phisum / dnum;
        _okPhi = true;
      }

    }


  private:

    /// The particle tracks.
    std::vector<FourMomentum> _particles;

    /// Cached values of \f$ \bar{\phi} \f$ and \f$ \bar{\eta} \f$.
    mutable double _phi, _eta;
    mutable bool _okPhi, _okEta;

    /// Cached values of the \f$ p_T \f$-weighted \f$ \bar{\phi} \f$ and \f$ \bar{\eta} \f$.
    mutable double _ptWeightedPhi, _ptWeightedEta;
    mutable bool _okPtWeightedPhi, _okPtWeightedEta;

    /// Cached value of the \f$ p_T \f$ sum.
    mutable double _totalPt;
    mutable bool _okTotalPt;
  };


  /// Typedef for a collection of Jet objects.
  typedef vector<Jet> Jets;

}

#endif
