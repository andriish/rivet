// -*- C++ -*-
#ifndef RIVET_KtJets_HH
#define RIVET_KtJets_HH
// Declaration of the KtJets class.

#include "Rivet/Projections/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "KtJet/KtJet.h"
#include "KtJet/KtEvent.h"


namespace Rivet {
  
  /// Project out all final-state particles in an event.
  class KtJets : public Projection {
    
  public:
    
    /// @name Standard constructors and destructors.
    //@{
    /// Default constructor. Must specify a FinalState projection which is
    //  assumed to live throughout the run.
    inline KtJets(FinalState& fsp)
      : type_(4), angle_(3), recom_(1), rparameter_(1.0), fsproj(&fsp) 
    { }

    /// Argument constructor.
    // Added so that same projection can be ran but with different parameters.
    // Must specify a FinalState projection which is
    // assumed to live throughout the run. 
    inline KtJets(FinalState& fsp, int type, int angle, int recom, double rparameter)
      : type_(type), angle_(angle), recom_(recom), rparameter_(rparameter), fsproj(&fsp)
    { }

    /// Copy constructor.
    inline KtJets(const KtJets& x)
      : Projection(x), type_(x.type_), angle_(x.angle_), 
        recom_(x.recom_), rparameter_(x.rparameter_), fsproj(x.fsproj)
    { }
    
    /// Destructor.
    virtual ~KtJets() { delete pktev_; };
    //@}

  protected:   

    /// Perform the projection on the Event: only to be called by 
    /// Event::applyProjection(Projection &).
    void project(const Event& e);

    /// This function defines a unique ordering between different 
    /// Projection objects of the same class. Should only be called from 
    /// operator<(const Projection&).
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

    /// Return the RivetInfo object of this Projection.
    virtual RivetInfo getInfo() const;

  private:
    
    // Internal ktevent, rebuilt every time an event is projected, but not otherwise.
    KtJet::KtEvent * pktev_;

    // Vector of all
    //vector<KtJet::KtLorentzVector> vecs_;

    int type_, angle_, recom_;
    double rparameter_;  // had trouble building with Ktfloat, used double instead

    /// The FinalState projection used by this projection.
    FinalState * fsproj;

  private:
    
    /// The assignment operator is private and must never be called.
    /// In fact, it shouldn't even be implemented.
    KtJets & operator=(const KtJets &);
  
  };

}

#endif /* RIVET_KtJets_HH */
