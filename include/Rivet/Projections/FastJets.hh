// -*- C++ -*-
#ifndef RIVET_FastJets_HH
#define RIVET_FastJets_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/JetAlg.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Jet.hh"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/SISConePlugin.hh"
//#include "fastjet/PxConePlugin.hh"
#include "fastjet/CDFJetCluPlugin.hh"
#include "fastjet/CDFMidPointPlugin.hh"
#ifdef HAVE_JADE
#include "fastjet/JadePlugin.hh"
#endif

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
    #ifdef HAVE_JADE
    enum JetAlg { KT, CAM, SISCONE, ANTIKT, PXCONE, CDFJETCLU, CDFMIDPOINT, JADE, DURHAM };
    #else
    enum JetAlg { KT, CAM, SISCONE, ANTIKT, PXCONE, CDFJETCLU, CDFMIDPOINT };
    #endif


    /// @name Constructors and destructors.
    //@{
    /// Default constructor uses \f$ k_\perp \f$ algorithm and E-scheme
    /// recombination.
    FastJets(const FinalState& fsp){
      setName("FastJets");
      addProjection(fsp, "FS");
      const double RPARAM = 1.0;
      _jdef = fastjet::JetDefinition(fastjet::kt_algorithm, RPARAM, fastjet::E_scheme);
    }


    /// "Wrapped" argument constructor using Rivet enums for most common
    /// jet alg choices (including some plugins). For the built-in algs,
    /// E-scheme recombination is used. For full control of
    /// FastJet built-in jet algs, use the native arg constructor.
    FastJets(const FinalState& fsp, JetAlg alg, double rparameter) {
      setName("FastJets");
      addProjection(fsp, "FS");
      if (alg == KT) {
        _jdef = fastjet::JetDefinition(fastjet::kt_algorithm, rparameter, fastjet::E_scheme);
      } else if (alg == CAM) {
        _jdef = fastjet::JetDefinition(fastjet::cambridge_algorithm, rparameter, fastjet::E_scheme);
      } else if (alg == ANTIKT) {
        _jdef = fastjet::JetDefinition(fastjet::antikt_algorithm, rparameter, fastjet::E_scheme);
      } else {
        if (alg == SISCONE) {
          const double overlapthreshold = 0.5;
          _plugin.reset(new fastjet::SISConePlugin(rparameter, overlapthreshold));
        } else if (alg == PXCONE) {
          throw Error("PxCone currently not supported, since FastJet doesn't install it by default");
          //_plugin.reset(new fastjet::PxConePlugin(rparameter));
        } else if (alg == CDFJETCLU) {
          _plugin.reset(new fastjet::CDFJetCluPlugin(rparameter));
        } else if (alg == CDFMIDPOINT) {
          _plugin.reset(new fastjet::CDFMidPointPlugin(rparameter));
        #ifdef HAVE_JADE
        } else if (alg == JADE) {
          _plugin.reset(new fastjet::JadePlugin("jade"));
        } else if (alg == DURHAM) {
          _plugin.reset(new fastjet::JadePlugin("durham"));
        #endif
        }
        _jdef = fastjet::JetDefinition(_plugin.get());
      }
    }
    
    
    /// Native argument constructor, using FastJet alg/scheme enums.
    FastJets(const FinalState& fsp, fastjet::JetAlgorithm type,
             fastjet::RecombinationScheme recom, double rparameter) {
      setName("FastJets");
      addProjection(fsp, "FS");
      _jdef = fastjet::JetDefinition(type, rparameter, recom);
    }


    /// Explicit copy constructor.
    FastJets(const FastJets& other) 
      : //_cseq(other._cseq),
        _jdef(other._jdef),
        _plugin(other._plugin),
        _yscales(other._yscales)
    {  
      setName("FastJets");
    }


    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new FastJets(*this);
    }

    //@}


  public:

    /// Reset the projection. Jet def etc are unchanged.
    void reset() { 
      _yscales.clear();
      _particles.clear();
      /// @todo Also clear the CSeq
      //_cseq.XXX;
    }

    /// Number of jets above the \f$ p_\perp \f$ cut.
    size_t numJets(double ptmin = 0.0) const {
      if (_cseq.get() != 0) {
        return _cseq->inclusive_jets(ptmin).size();
      } else {
        return 0;
      }        
    }

    /// Number of jets.
    size_t size() const {
      return numJets();
    }

    /// Get the jets (unordered).
    Jets jets(double ptmin = 0.0) const {
      return _pseudojetsToJets(pseudoJets(ptmin));
    }
    
    /// Get the jets, ordered by \f$ p_T \f$.
    Jets jetsByPt(double ptmin = 0.0) const {
      return _pseudojetsToJets(pseudoJetsByPt(ptmin));
    }

    /// Get the jets, ordered by \f$ E \f$.
    Jets jetsByE(double ptmin = 0.0) const {
      return _pseudojetsToJets(pseudoJetsByE(ptmin));
    }

    /// Get the jets, ordered by rapidity.
    Jets jetsByRapidity(double ptmin = 0.0) const {
      return _pseudojetsToJets(pseudoJetsByRapidity(ptmin));
    }

    /// Get the pseudo jets (unordered).
    PseudoJets pseudoJets(double ptmin = 0.0) const {
      if (_cseq.get() != 0) {
        return _cseq->inclusive_jets(ptmin);
      } else {
        return PseudoJets();
      }
    }

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
