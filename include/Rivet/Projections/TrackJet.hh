// -*- C++ -*-
#ifndef RIVET_TrackJet_HH
#define RIVET_TrackJet_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"


// *** Only count hits where SUMPT is different from 0! 


namespace Rivet {

  
  /// Project out jets of charged tracks a la CDF.
  class TrackJet : public Projection {
  public:

    /// @name Standard constructors and destructors.
    //@{
    /// Constructor. The specified FinalState projection is assumed to live
    /// throughout the run and should be used to specify the max and min \f$
    /// \eta \f$ values and the min \f$ p_T \f$ (in GeV).
    inline TrackJet(FinalState& fsp)
      : _fsproj(fsp)
    { 
      addProjection(fsp);
    }

    /// Typedef for the tracks (a list so that elements can be consistently removed
    typedef list<FourMomentum> Tracks;


  public:    

    /// Inner class used to represent a jet as defined by this algorithm.
    class Jet {
    public:
      /// Constructor.
      inline Jet() { 
        clear();
      }

      /// Define a Jet::iterator via a typedef.
      typedef vector<FourMomentum>::iterator iterator;

      /// Define a Jet::const_iterator via a typedef.
      typedef vector<FourMomentum>::const_iterator const_iterator;

      /// Get a begin iterator over the particles/tracks in this jet.
      inline iterator begin() {
        return _particles.begin();
      }

      /// Get an end iterator over the particles/tracks in this jet.
      inline iterator end() {
        return _particles.end();
      }

      /// Get a const begin iterator over the particles/tracks in this jet.
      inline const_iterator begin() const {
        return _particles.begin();
      }

      /// Get a const end iterator over the particles/tracks in this jet.
      inline const_iterator end() const {
        return _particles.end();
      }

      /// Get the particles (tracks) in this jet.
      inline vector<FourMomentum>& getParticles() {
        return _particles;
      }

      /// Get the particles (tracks) in this jet (const version).
      inline const vector<FourMomentum>& getParticles() const {
        return _particles;
      }

      /// Set the particles/tracks collection.
      inline Jet setParticles(vector<FourMomentum> particles) {
        _particles = particles;
        _resetCaches();
        return *this;
      }

      /// Add a particle/track to this jet.
      inline Jet addParticle(FourMomentum particle) {
        _particles.push_back(particle);
        _resetCaches();
        return *this;
      }

      /// Reset this jet as empty.
      inline Jet clear() {
        _particles.clear();
        _resetCaches();
        return *this;
      }

      /// Get the average \f$ \phi \f$ for this jet, with the average weighted
      /// by the \f$ p_T \f$ values of the constituent tracks. (caches)
      inline double getPtWeightedPhi() const {
        if (_ptWeightedPhi < 0) {
          double ptwphi(0.0), ptsum(0.0);
          for (const_iterator p = this->begin(); p != this->end(); ++p) {
            double pt = pT(*p);
            ptwphi += pt * p->vector3().azimuthalAngle();
            ptsum += pt;
          }
          _totalPt = ptsum;
          _ptWeightedPhi = ptwphi / getPtSum();
        }
        return _ptWeightedPhi;
      }

      /// Get the sum of the \f$ p_T \f$ values of the constituent tracks. (caches)
      inline double getPtSum() const {
        if (_totalPt < 0) {
          double ptsum(0.0);
          for (const_iterator p = this->begin(); p != this->end(); ++p) {
            ptsum += pT(*p);
          }
          _totalPt = ptsum;
        }
        return _totalPt;
      }

      /// Get the number of particles/tracks in this jet.
      inline size_t getNumParticles() const {
        return _particles.size();
      }


    private:

      /// Clear the internal cached values.
      inline Jet _resetCaches() {
        _totalPt = -1.0;
        _ptWeightedPhi = -1.0;
        return *this;
      }

      /// The particle tracks.
      vector<FourMomentum> _particles;

      /// Cached value of the \f$ p_T \f$-weighted \f$ \bar{\phi} \f$
      mutable double _ptWeightedPhi;

      /// Cached value of the \f$ p_T \f$ sum.
      mutable double _totalPt;
    };


    /// Typedef for a collection of Jet objects.
    typedef vector<Jet> Jets;


  public:
    /// Return the name of the projection
    inline string getName() const {
      return "TrackJet";
    }

    /// Get the computed jets.
    inline Jets& getJets() {
      return _jets;
    }

    /// Get the computed jets (const version).
    inline const Jets& getJets() const {
      return _jets;
    }

  protected:   

    /// Perform the projection on the Event. The collection of jets that results
    /// will be sorted in order of decreasing \f$ p_T \f$.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;  

  private:
    
    /// The FinalState projection used by this projection.
    FinalState _fsproj;

    /// The computed jets
    Jets _jets;
  
  };


  /// Hide helper functions in an anonymous namespace
  namespace {    
    /// For sorting Jet objects by pT.
    inline bool compareJetsByPt(const TrackJet::Jet& first, const TrackJet::Jet& second) {
      return first.getPtSum() > second.getPtSum();
    }

    /// For sorting Lorentz four-vectors by pT.
    inline bool compareVecsByPt(const FourMomentum& first, const FourMomentum& second) {
      return pT2(first) > pT2(second);
    }
  }
  
 
}

#endif
