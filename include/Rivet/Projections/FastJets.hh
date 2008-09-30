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
    /// Return the number of jets above the pt cut 
    size_t getNumJets(double ptmin = 0.0) const {
      return _cseq.inclusive_jets(ptmin).size();
    }

    /// Get the jets (unordered).
    Jets getJets(double ptmin = 0.0) const {
      return _pseudojetsToJets(getPseudoJets(ptmin));
    }
    
    /// Get the jets, ordered by \f$ p_T \f$.
    Jets getJetsByPt(double ptmin = 0.0) const {
      return _pseudojetsToJets(getPseudoJetsByPt(ptmin));
    }

    /// Get the jets, ordered by \f$ E \f$.
    Jets getJetsByE(double ptmin = 0.0) const {
      return _pseudojetsToJets(getPseudoJetsByE(ptmin));
    }

    /// Get the jets, ordered by rapidity.
    Jets getJetsByRapidity(double ptmin = 0.0) const {
      return _pseudojetsToJets(getPseudoJetsByRapidity(ptmin));
    }

    /// Get the pseudo jets (unordered).
    PseudoJets getPseudoJets(double ptmin = 0.0) const {
      return _cseq.inclusive_jets(ptmin);
    }

    /// Get the pseudo jets, ordered by \f$ p_T \f$.
    PseudoJets getPseudoJetsByPt(double ptmin = 0.0) const {
      return sorted_by_pt(_cseq.inclusive_jets(ptmin));
    }

    /// Get the pseudo jets, ordered by \f$ E \f$.
    PseudoJets getPseudoJetsByE(double ptmin = 0.0) const {
      return sorted_by_E(_cseq.inclusive_jets(ptmin));
    }

    /// Get the pseudo jets, ordered by rapidity.
    PseudoJets getPseudoJetsByRapidity(double ptmin = 0.0) const {
      return sorted_by_rapidity(_cseq.inclusive_jets(ptmin));
    }

    /// Return the cluster sequence (FastJet-specific).
    const fastjet::ClusterSequence& getClusterSeq() const {
      return _cseq;
    }

    /// Return the jet definition (FastJet-specific).
    const fastjet::JetDefinition& getJetDef() const {
      return _jdef;
    }

    /// Get the subjet splitting variables for the given jet.
    vector<double> getYSubJet(const fastjet::PseudoJet& jet) const;

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


  private:

    /// Jet definition
    fastjet::JetDefinition _jdef;

    /// Cluster sequence
    fastjet::ClusterSequence _cseq;

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
