// -*- C++ -*-
#ifndef RIVET_TrackJet_HH
#define RIVET_TrackJet_HH

#include "Rivet/Projections/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/RivetCLHEP.hh"


// *** Only count hits where SUMPT is different from 0! 


namespace Rivet {

  namespace {
    /// Calculate \f$ p_T^2 \f$ for the provided LorentzVector.
    inline double pT2(const LorentzVector& lv) {
      return lv.x()*lv.x() + lv.y()*lv.y();
    }
    
    /// Calculate \f$ p_T \f$ for the provided LorentzVector.
    inline double pT(const LorentzVector& lv) {
      return sqrt(pT2(lv));
    }
  }

  
  /// Project out all final-state particles in an event.
  class TrackJet : public Projection {

  public:    

    /// Inner class used to represent a jet as defined by this algorithm.
    class Jet {
    public:
      /// Constructor.
      inline Jet() { 
        clear();
      }

      /// Define a Jet::iterator via a typedef.
      typedef vector<LorentzVector>::iterator iterator;

      /// Define a Jet::const_iterator via a typedef.
      typedef vector<LorentzVector>::const_iterator const_iterator;

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
      inline vector<LorentzVector>& getParticles() {
        return _particles;
      }

      /// Get the particles (tracks) in this jet (const version).
      inline const vector<LorentzVector>& getParticles() const {
        return _particles;
      }

      /// Set the particles/tracks collection.
      inline Jet setParticles(vector<LorentzVector> particles) {
        _particles = particles;
        _resetCaches();
        return *this;
      }

      /// Add a particle/track to this jet.
      inline Jet addParticle(LorentzVector particle) {
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
            ptwphi += pt * p->phi();
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
      vector<LorentzVector> _particles;

      /// Cached value of the \f$ p_T \f$-weighted \f$ \bar{\phi} \f$
      mutable double _ptWeightedPhi;

      /// Cached value of the \f$ p_T \f$ sum.
      mutable double _totalPt;
    };


    /// Typedef for the tracks (a list so that elements can be consistently removed
    typedef list<LorentzVector> Tracks;


    /// Typedef for a collection of Jet objects.
    typedef vector<Jet> Jets;

  public:

    /// @name Standard constructors and destructors.
    //@{
    /// Constructor. The specified FinalState projection is assumed to live
    /// throughout the run and should be used to specify the max and min \f$
    /// \eta \f$ values and the min \f$ p_T \f$ (in GeV).
    inline TrackJet(FinalState& fsp)
      : _fsproj(&fsp)
    { }

//     /// Argument constructor. 
//     /// The specified FinalState projection is assumed to live throughout the run. 
//     inline TrackJet(FinalState& fsp, const double etaMin, const double etaMax, const double ptMin)
//       : _fsproj(&fsp)
//     { }
    //@}

  public:
    /// Return the name of the projection
    inline string name() const {
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
    FinalState* _fsproj;

//     /// The minimum value of \f$ \eta \f$.
//     const double _etaMin;

//     /// The maximum value of \f$ \eta \f$.    
//     const double _etaMax;

//     /// The minimum value of \f$ p_T \f$ in GeV.
//     const double _ptMin;

    /// The computed jets
    Jets _jets;

  private:
    
    /// Hiding the assignment operator.
    TrackJet& operator=(const TrackJet&);
  
  };


  /// Hide helper functions in an anonymous namespace
  namespace {
    /// Sorts Lorentz four-vectors by pT.
    inline bool compareVecsByPt(const LorentzVector& first, const LorentzVector& second) {
      return pT2(first) > pT2(second);
    }
    
    /// Sorts Jet objects by pT.
    inline bool compareJetsByPt(const TrackJet::Jet& first, const TrackJet::Jet& second) {
      return first.getPtSum() > second.getPtSum();
    }
  }
  

  
}

#endif
