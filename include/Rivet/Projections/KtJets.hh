// -*- C++ -*-
#ifndef RIVET_KtJets_HH
#define RIVET_KtJets_HH

#include "Rivet/Projections/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "KtJet/KtJet.h"
#include "KtJet/KtEvent.h"


namespace Rivet {
  
  /// Project out jets based on configurable kT algorithm.
  class KtJets : public Projection {
    
  public:
    
    /// @name Standard constructors and destructors.
    //@{
    /// Default constructor. Must specify a FinalState projection which is
    //  assumed to live throughout the run.
    inline KtJets(FinalState& fsp)
      : _pktev(0), _type(4), _angle(2), _recom(1), 
        _rparameter(1.0),	_fsproj(&fsp) 
    { 
      addProjection(fsp);
    }

    /// Argument constructor. Allows the to be run with different parameters.
    /// Must specify a FinalState projection which is assumed to live throughout the run. 
    inline KtJets(FinalState& fsp, int type, int angle, int recom, double rparameter)
      : _pktev(0), _type(type), _angle(angle), _recom(recom),
        _rparameter(rparameter), _fsproj(&fsp)
    { 
      addProjection(fsp);
    }
    
    /// Destructor.
    inline virtual ~KtJets() { 
      if (_pktev) delete _pktev; 
    }
    //@}

  public:
    /// Return the name of the projection
    inline string getName() const {
      return "KtJets";
    }
    
  protected:   

    /// Perform the projection on the Event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;  

  public:
    
    /// @name Access the projected NJets.
    //@ {
    inline int getNJets() const { return _pktev->getNJets(); }
    inline int getNConstituents() const { return _pktev->getNConstituents(); }
    inline vector<KtJet::KtLorentzVector> copyConstituents() const { return _pktev->copyConstituents(); }
    inline double getETot() const { return _pktev->getETot(); } // had trouble building with Ktfloat, used double instead
    inline int getType() const { return _pktev->getType(); }
    inline int getAngle() const { return _pktev->getAngle(); }
    inline int getRecom() const { return _pktev->getRecom(); }
    inline bool isInclusive() const { return _pktev->isInclusive(); }

    inline vector<KtJet::KtLorentzVector> getJets() const { return _pktev->getJets(); }
    inline vector<KtJet::KtLorentzVector> getJetsE() const { return _pktev->getJetsE(); }
    inline vector<KtJet::KtLorentzVector> getJetsEt() const { return _pktev->getJetsEt(); }
    inline vector<KtJet::KtLorentzVector> getJetsPt() const { return _pktev->getJetsPt(); }
    inline vector<KtJet::KtLorentzVector> getJetsRapidity() const { return _pktev->getJetsRapidity(); }
    inline vector<KtJet::KtLorentzVector> getJetsEta() const { return _pktev->getJetsEta(); }
    //@}

    /// Get the subjet splitting variables for the given jet.
    vector<double> getYSubJet(const KtJet::KtLorentzVector& jet) const; 

  private:
    
    /// Internal KtEvent, rebuilt every time an event is projected, but not otherwise.
    KtJet::KtEvent* _pktev;

    int _type, _angle, _recom;
    double _rparameter;  // had trouble building with Ktfloat, used double instead

    /// The FinalState projection used by this projection.
    FinalState* _fsproj;

    /// Map of vectors of y scales. This is mutable so we can use caching/lazy evaluation.
    mutable map<int, vector<double> > _yscales;

  private:
    
  };

}

#endif
