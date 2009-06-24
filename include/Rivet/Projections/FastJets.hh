// -*- C++ -*-
#ifndef RIVET_FastJets_HH
#define RIVET_FastJets_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/JetAlg.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Jet.hh"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

namespace Rivet {


  /// Make a 3-momentum vector from a FastJet pseudo-jet
  inline Vector3 momentum3(const fastjet::PseudoJet& pj) {
    return Vector3(pj.px(), pj.py(), pj.pz());
  }

  /// Make a 4-momentum vector from a FastJet pseudo-jet
  inline FourMomentum momentum(const fastjet::PseudoJet& pj) {
    return FourMomentum(pj.E(), pj.px(), pj.py(), pj.pz());
  }
  

  /// Typedef for a collection of PseudoJets.
  typedef vector<fastjet::PseudoJet> PseudoJets;
  

  /// Project out jets found using fastJet package.
  class FastJets : public JetAlg {

  public:
    /// Wrapper enum for selected Fastjet jet algorithms.
    enum JetAlgName { KT, CAM, SISCONE, ANTIKT, PXCONE, 
                  CDFJETCLU, CDFMIDPOINT, D0ILCONE,
                  JADE, DURHAM, TRACKJET };


    /// @name Constructors etc.
    //@{

    /// "Wrapped" argument constructor using Rivet enums for most common
    /// jet alg choices (including some plugins). For the built-in algs,
    /// E-scheme recombination is used. For full control of
    /// FastJet built-in jet algs, use the native arg constructor.
    FastJets(const FinalState& fsp, JetAlgName alg, 
             double rparameter, double pTmin=0.0, double seed_threshold=1.0);

    /// Native argument constructor, using FastJet alg/scheme enums.
    FastJets(const FinalState& fsp, fastjet::JetAlgorithm type,
             fastjet::RecombinationScheme recom, double rparameter);

    /// Explicitly pass in an externally-constructed plugin
    FastJets(const FinalState& fsp, const fastjet::JetDefinition::Plugin& plugin);

    /// Explicit copy constructor.
    FastJets(const FastJets& other);

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new FastJets(*this);
    }

    //@}


  public:

    /// Reset the projection. Jet def etc are unchanged.
    void reset();

    /// Number of jets above the \f$ p_\perp \f$ cut.
    size_t numJets(double ptmin = 0.0) const;

    /// Number of jets.
    size_t size() const {
      return numJets();
    }

    /// Get the jets (unordered).
    Jets jets(double ptmin = 0.0) const;
    
    /// Get the jets, ordered by \f$ p_T \f$.
    Jets jetsByPt(double ptmin = 0.0) const;

    /// Get the jets, ordered by \f$ E \f$.
    Jets jetsByE(double ptmin = 0.0) const;

    /// Get the jets, ordered by rapidity.
    Jets jetsByRapidity(double ptmin = 0.0) const;

    /// Get the pseudo jets (unordered).
    PseudoJets pseudoJets(double ptmin = 0.0) const;

    /// Get the pseudo jets, ordered by \f$ p_T \f$.
    PseudoJets pseudoJetsByPt(double ptmin = 0.0) const {
      return sorted_by_pt(pseudoJets(ptmin));
    }

    /// Get the pseudo jets, ordered by \f$ E \f$.
    PseudoJets pseudoJetsByE(double ptmin = 0.0) const {
      return sorted_by_E(pseudoJets(ptmin));
    }

    /// Get the pseudo jets, ordered by rapidity.
    PseudoJets pseudoJetsByRapidity(double ptmin = 0.0) const {
      return sorted_by_rapidity(pseudoJets(ptmin));
    }

    /// Return the cluster sequence (FastJet-specific).
    const fastjet::ClusterSequence* clusterSeq() const {
      return _cseq.get();
    }

    /// Return the jet definition (FastJet-specific).
    const fastjet::JetDefinition& jetDef() const {
      return _jdef;
    }

    /// Get the subjet splitting variables for the given jet.
    vector<double> ySubJet(const fastjet::PseudoJet& jet) const;

    /// @brief Split a jet a la PRL100,242001(2008).
    /// Based on code from G.Salam, A.Davison.
    fastjet::PseudoJet splitJet(fastjet::PseudoJet jet, double& last_R) const;

    /// @brief Filter a jet a la PRL100,242001(2008).
    /// Based on code from G.Salam, A.Davison.
    fastjet::PseudoJet filterJet(fastjet::PseudoJet jet, double& stingy_R, const double def_R) const;

  private:

    Jets _pseudojetsToJets(const PseudoJets& pjets) const;
      
  protected:   

    /// Perform the projection on the Event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;  

  public:

    /// Do the calculation locally (no caching).
    void calc(const ParticleVector& ps);


  private:

    /// Jet definition
    fastjet::JetDefinition _jdef;

    /// Cluster sequence
    shared_ptr<fastjet::ClusterSequence> _cseq;

    /// FastJet external plugin
    shared_ptr<fastjet::JetDefinition::Plugin> _plugin; 
    
    /// Map of vectors of y scales. This is mutable so we can use caching/lazy evaluation.
    mutable map<int, vector<double> > _yscales;
  
    /// set of particles sorted by their PT2
    //set<Particle, ParticleBase::byPTAscending> _particles;
    map<int, Particle> _particles;
    
  };

}

#endif
