// -*- C++ -*-
#ifndef RIVET_KtJets_HH
#define RIVET_KtJets_HH
// Declaration of the KtJets class.

#include "Rivet/Projections/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "KtJet/KtJet.h"

namespace Rivet {
  
  /**
   * Project out all final-state particles in an event.
   */
  class KtJets : public Projection {
    
  public:
    
    /// @name Standard constructors and destructors.
    //@{
    /// Default constructor.
    inline KtJets();
    
    /// Copy constructor.
    inline KtJets(const KtJets &);
    
    /// Destructor.
    virtual ~KtJets();
    //@}

protected:   

    /// Perform the projection on the Event: only to be called by 
    /// Event::addProjection(Projection &).
    void project(const Event & e);

    /// This function defines a unique ordering between different 
    /// Projection objects of the same class. Should only be called from 
    /// operator<(const Projection &).
    int compare(const Projection & p) const;  

  public:
    
    /// @name Access the projected NJets.
    //@ {
    /// NJets
    const unsigned int NJets() const;
    void getJetsE() ;
 
  private:
    
    /// Total NJets
    unsigned int NJets_;
    std::vector<KtJet::KtLorentzVector> jetVecs_;

    /// The FinalState projection used by this projection
    FinalState fsproj;
 
  private:
    
    /**
     * The assignment operator is private and must never be called.
     * In fact, it shouldn't even be implemented.
     */
    KtJets & operator=(const KtJets &);
    
  };

}

#include "Rivet/Projections/KtJets.icc"

#endif /* RIVET_KtJets_HH */
