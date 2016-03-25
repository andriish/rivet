// -*- C++ -*-
#ifndef RIVET_SmearedJets_HH
#define RIVET_SmearedJets_HH

#include "Rivet/Jet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/JetAlg.hh"
#include "Rivet/Tools/SmearingFunctions.hh"
#include <functional>

namespace Rivet {


  /// Wrapper projection for smearing {@link Jet}s with detector resolutions and efficiencies
  /// @todo Chain constructors?
  class SmearedJets : public JetAlg {
  public:

    /// @name Constructors etc.
    //@{

    /// @brief Constructor with efficiency and smearing function args
    /// The jet reconstruction efficiency is mandatory; the smearing and tagging functions are optional
    template <typename J2JFN>
    SmearedJets(const JetAlg& ja,
                const J2JFN& jetSmearFn)
      : _jetEffFn(JET_EFF_ONE), _bTagEffFn(JET_BTAG_PERFECT), _cTagEffFn(JET_CTAG_PERFECT), _jetSmearFn(jetSmearFn)
    {
      setName("SmearedJets");
      addProjection(ja, "TruthJets");
    }


    /// @brief Constructor with efficiency and smearing function args
    /// The jet reconstruction efficiency is mandatory; the smearing and tagging functions are optional
    template <typename J2JFN, typename J2DFN>
    SmearedJets(const JetAlg& ja,
                const J2JFN& jetSmearFn,
                const J2DFN& bTagEffFn)
      : _jetEffFn(JET_EFF_ONE), _bTagEffFn(bTagEffFn), _cTagEffFn(JET_CTAG_PERFECT), _jetSmearFn(jetSmearFn)
    {
      setName("SmearedJets");
      addProjection(ja, "TruthJets");
    }


    /// @brief Constructor with efficiency and smearing function args
    /// The jet reconstruction efficiency is mandatory; the smearing and tagging functions are optional
    template <typename J2JFN, typename J2DFN>
    SmearedJets(const JetAlg& ja,
                const J2JFN& jetSmearFn,
                const J2DFN& bTagEffFn,
                const J2DFN& cTagEffFn)
      : _jetEffFn(JET_EFF_ONE), _bTagEffFn(bTagEffFn), _cTagEffFn(cTagEffFn), _jetSmearFn(jetSmearFn)
    {
      setName("SmearedJets");
      addProjection(ja, "TruthJets");
    }


    /// @brief Constructor with efficiency and smearing function args
    /// The jet reconstruction efficiency is mandatory; the smearing and tagging functions are optional
    template <typename J2JFN, typename J2DFN>
    SmearedJets(const JetAlg& ja,
                const J2JFN& jetSmearFn,
                const J2DFN& bTagEffFn,
                const J2DFN& cTagEffFn,
                const J2DFN& jetEffFn)
      : _jetEffFn(jetEffFn), _bTagEffFn(bTagEffFn), _cTagEffFn(cTagEffFn), _jetSmearFn(jetSmearFn)
    {
      setName("SmearedJets");
      addProjection(ja, "TruthJets");
    }


    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new SmearedJets(*this);
    }

    //@}


    /// Compare to another SmearedJets
    int compare(const Projection& p) const {
      const SmearedJets& other = dynamic_cast<const SmearedJets&>(p);
      if (_mkhash(_jetEffFn) == 0) return UNDEFINED;
      if (_mkhash(_bTagEffFn) == 0) return UNDEFINED;
      if (_mkhash(_cTagEffFn) == 0) return UNDEFINED;
      if (_mkhash(_jetSmearFn) == 0) return UNDEFINED;
      return
        cmp(_mkhash(_jetEffFn), _mkhash(other._jetEffFn)) || cmp(_mkhash(_jetSmearFn), _mkhash(other._jetSmearFn)) ||
        cmp(_mkhash(_bTagEffFn), _mkhash(other._bTagEffFn)) || cmp(_mkhash(_cTagEffFn), _mkhash(other._cTagEffFn));
    }


    /// Perform the jet finding & smearing calculation
    void project(const Event& e) {
      // Copying and filtering
      const Jets& truthjets = applyProjection<JetAlg>(e, "TruthJets").jetsByPt();
      _recojets.clear(); _recojets.reserve(truthjets.size());
      foreach (const Jet& j, truthjets) {
        const double jeff = (_jetEffFn) ? _jetEffFn(j) : 1;
        MSG_DEBUG("Efficiency of jet " << j.mom() << " = " << 100*jeff << "%");
        if (jeff == 0) continue; //< no need to roll expensive dice
        if (jeff == 1 || jeff < rand01()) {
          _recojets.push_back(_jetSmearFn ? _jetSmearFn(j) : j); //< smearing
        }
      }
      // Tagging efficiencies
      foreach (Jet& j, _recojets) {
        const double beff = (_bTagEffFn) ? _bTagEffFn(j) : 1;
        const bool btag = beff == 1 || (beff != 0 && beff < rand01());
        // Remove b-tags if needed, and add a dummy one if needed
        if (!btag && j.bTagged()) j.tags().erase(std::remove_if(j.tags().begin(), j.tags().end(), hasBottom), j.tags().end());
        if (btag && !j.bTagged()) j.tags().push_back(Particle(PID::BQUARK, j.mom()));
        const double ceff = (_cTagEffFn) ? _cTagEffFn(j) : 1;
        const bool ctag = ceff == 1 || (ceff != 0 && ceff < rand01());
        // Remove c-tags if needed, and add a dummy one if needed
        if (!ctag && j.cTagged()) j.tags().erase(std::remove_if(j.tags().begin(), j.tags().end(), hasCharm), j.tags().end());
        if (ctag && !j.cTagged()) j.tags().push_back(Particle(PID::CQUARK, j.mom()));
      }
    }


    /// Return the full jet list for the JetAlg methods to use
    Jets _jets() const { return _recojets; }


    /// Reset the projection. Smearing functions will be unchanged.
    void reset() { _recojets.clear(); }


  private:

    size_t _mkhash(const std::function<double(const Jet&)>& fn) const {
      return reinterpret_cast<size_t>(fn.target<double(*)(const Jet&)>());
    }

    size_t _mkhash(const std::function<Jet(const Jet&)>& fn) const {
      return reinterpret_cast<size_t>(fn.target<Jet(*)(const Jet&)>());
    }


    Jets _recojets;
    std::function<double(const Jet&)> _jetEffFn, _bTagEffFn, _cTagEffFn;
    std::function<Jet(const Jet&)> _jetSmearFn;

  };


}

#endif
