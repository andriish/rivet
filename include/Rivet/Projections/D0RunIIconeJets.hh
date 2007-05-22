// -*- C++ -*-
#ifndef RIVET_D0RunIIconeJets_HH
#define RIVET_D0RunIIconeJets_HH
// Declaration of the D0RunIIconeJets class.

#include "Rivet/Rivet.hh"
#include "Rivet/Projections/Projection.hh"
#include "Rivet/Projections/FinalState.hh"

#include "Rivet/Tools/D0RunIIcone/HepEntity.h"
#include "Rivet/Tools/D0RunIIcone/energycluster/ILConeAlgorithm.hpp"


namespace Rivet {
  
  /// Project out all final-state particles in an event.
  class D0RunIIconeJets : public Projection {
    
  public:
    
    /// @name Standard constructors and destructors.
    //@{
    /// Default constructor. Must specify a FinalState projection which is
    //  assumed to live throughout the run.
    inline D0RunIIconeJets(FinalState& fsp)
      : cone_radius(0.7), min_jet_Et(0.), split_ratio(0.5), fsproj(&fsp) 
    { 
      addProjection(fsp);

      // The parameters are supposed to be set as used by D0 in RunII
      far_ratio = 0.5;
      Et_min_ratio = 0.5;
      kill_duplicate = true;
      duplicate_dR = 0.005; 
      duplicate_dPT = 0.01; 
      search_factor = 1.0; 
      pT_min_leading_protojet = 0.0; 
      pT_min_second_protojet = 0.0;
      merge_max = 10000;
      pT_min_nomerge = 0.0;
      
      algo = new ILConeAlgorithm<HepEntity>(cone_radius, min_jet_Et, split_ratio,
		      far_ratio, Et_min_ratio, kill_duplicate, duplicate_dR, 
		      duplicate_dPT, search_factor, pT_min_leading_protojet, 
		      pT_min_second_protojet, merge_max, pT_min_nomerge);
      
      jets = new std::list<HepEntity>;
    }

    /// Argument constructor.
    // Added so that same projection can be ran but with different parameters.
    // Must specify a FinalState projection which is
    // assumed to live throughout the run. 
    inline D0RunIIconeJets(FinalState& fsp, float R, float Etmin, float split)
      : cone_radius(R), min_jet_Et(Etmin), split_ratio(split), fsproj(&fsp)
    { 
      algo = new ILConeAlgorithm<HepEntity>(cone_radius, min_jet_Et, split_ratio,
					      far_ratio, Et_min_ratio, kill_duplicate, duplicate_dR, 
					      duplicate_dPT, search_factor, pT_min_leading_protojet, 
					      pT_min_second_protojet, merge_max, pT_min_nomerge);
      
      jets = new std::list<HepEntity>;
    }

    /// Copy constructor.
    inline D0RunIIconeJets(const D0RunIIconeJets& x)
      : Projection(x), cone_radius(x.cone_radius), min_jet_Et(x.min_jet_Et), 
	split_ratio(x.split_ratio), fsproj(x.fsproj)
    { 
      algo = new ILConeAlgorithm<HepEntity>(cone_radius, min_jet_Et, split_ratio,
	     far_ratio, Et_min_ratio, kill_duplicate, duplicate_dR, 
	     duplicate_dPT, search_factor, pT_min_leading_protojet, 
	     pT_min_second_protojet, merge_max, pT_min_nomerge);
      
      jets->clear();
      for (std::list<HepEntity>::iterator it=x.jets->begin(); it!=x.jets->end(); ++it)
	jets->push_back(*it);
  
    }
    
    /// Destructor.
    virtual ~D0RunIIconeJets() { 
      delete particlelist; 
      delete algo;
      delete jets;
    };
    //@}

  public:
    /// Return the name of the projection
    inline string getName() const {
      return "D0RunIIconeJets";
    }

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
    /*
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
    */
    //@}

    inline int getNJets() const { return jets->size(); }



    /// Return the RivetInfo object of this Projection.
    // virtual RivetInfo getInfo() const;

    std::list<HepEntity> * jets;
    
    
  private:
    
    /// The assignment operator is private and must never be called.
    /// In fact, it shouldn't even be implemented.
    D0RunIIconeJets & operator=(const D0RunIIconeJets &);
  
    // Internal ktevent, rebuilt every time an event is projected, but not otherwise.
    //KtJet::KtEvent * pktev_;

    // Vector of all
    //vector<KtJet::KtLorentzVector> vecs_;
    std::list<const HepEntity*> *particlelist;
    
    //initialize D0RunII cone algorithm
    float cone_radius;
    float min_jet_Et;
    float split_ratio;
    
    float far_ratio;
    float Et_min_ratio;
    bool kill_duplicate;
    float duplicate_dR; 
    float duplicate_dPT; 
    float search_factor; 
    float pT_min_leading_protojet; 
    float pT_min_second_protojet;
    int merge_max; 
    float pT_min_nomerge;

    /// The FinalState projection used by this projection.
    FinalState * fsproj;

    ILConeAlgorithm<HepEntity> * algo;

  };

}

#endif /* RIVET_D0RunIIconeJets_HH */
