// -*- C++ -*-
#ifndef RIVET_KtJets_HH
#define RIVET_KtJets_HH
// Declaration of the KtJets class.

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
      : pktev_(0),type_(4), angle_(3), recom_(1), rparameter_(1.0),
	fsproj(&fsp) 
    { 
      addProjection(fsp);
    }

    /// Argument constructor. Allows the to be run with different parameters.
    /// Must specify a FinalState projection which is assumed to live throughout the run. 
    inline KtJets(FinalState& fsp, int type, int angle, int recom, double rparameter)
      : pktev_(0), type_(type), angle_(angle), recom_(recom),
	rparameter_(rparameter), fsproj(&fsp)
    { 
      addProjection(fsp);
    }
    
    /// Destructor.
    inline virtual ~KtJets() { if ( pktev_ ) delete pktev_; };
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
    inline int getNJets() const { return pktev_->getNJets(); }
    inline int getNConstituents() const { return pktev_->getNConstituents(); }
    inline vector<KtJet::KtLorentzVector> copyConstituents() const { return pktev_->copyConstituents(); }
    inline double getETot() const { return pktev_->getETot(); } // had trouble building with Ktfloat, used double instead
    inline int getType() const { return pktev_->getType(); }
    inline int getAngle() const { return pktev_->getAngle(); }
    inline int getRecom() const { return pktev_->getRecom(); }
    inline bool isInclusive() const { return pktev_->isInclusive(); }

    inline vector<KtJet::KtLorentzVector> getJets() const { return pktev_->getJets(); }
    inline vector<KtJet::KtLorentzVector> getJetsE() const { return pktev_->getJetsE(); }
    inline vector<KtJet::KtLorentzVector> getJetsEt() const { return pktev_->getJetsEt(); }
    inline vector<KtJet::KtLorentzVector> getJetsPt() const { return pktev_->getJetsPt(); }
    inline vector<KtJet::KtLorentzVector> getJetsRapidity() const { return pktev_->getJetsRapidity(); }
    inline vector<KtJet::KtLorentzVector> getJetsEta() const { return pktev_->getJetsEta(); }
    //@}

    /// Get the subjet splitting variables for the given jet.
    vector<double> getYSubJet(const KtJet::KtLorentzVector &jet) const; 

  private:
    
    /// Internal KtEvent, rebuilt every time an event is projected, but not otherwise.
    KtJet::KtEvent* pktev_;

    int type_, angle_, recom_;
    double rparameter_;  // had trouble building with Ktfloat, used double instead

    /// The FinalState projection used by this projection.
    FinalState* fsproj;

    /// Map of vectors of y scales. This is mutable so we can use caching/lazy evaluation.
    mutable map<int,vector<double> > yscales_;

  private:
    
    /// Hiding the assignment operator.
    KtJets& operator=(const KtJets &);
  
  };

}

#endif /* RIVET_KtJets_HH */
