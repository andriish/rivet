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
  class SmearedJets : public JetAlg {
  public:

    /// @name Constructors etc.
    //@{

    /// @brief Constructor with efficiency and smearing function args
    /// The jet reconstruction efficiency is mandatory; the smearing and tagging functions are optional
    template <typename J2DFN, typename J2JFN>
    SmearedJets(const JetAlg& ja,
                const J2JFN& jetSmearFn=JET_SMEAR_IDENTITY,
                const J2DFN& bTagEffFn=JET_BTAG_PERFECT, const J2DFN& cTagEffFn=JET_CTAG_PERFECT,
                const J2DFN& jetEffFn=JET_EFF_ONE)
      : _jetEffFnHash(reinterpret_cast<size_t>(jetEffFn)),
        _bTagEffFnHash(reinterpret_cast<size_t>(bTagEffFn)),
        _cTagEffFnHash(reinterpret_cast<size_t>(cTagEffFn)),
        _jetSmearFnHash(reinterpret_cast<size_t>(jetSmearFn)),
        _jetEffFn(jetEffFn), _bTagEffFn(bTagEffFn), _cTagEffFn(cTagEffFn), _jetSmearFn(jetSmearFn)
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
      return
        cmp(_jetEffFnHash, other._jetEffFnHash) || cmp(_jetSmearFnHash, other._jetSmearFnHash) ||
        cmp(_bTagEffFnHash, other._bTagEffFnHash) || cmp(_cTagEffFnHash, other._cTagEffFnHash);
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

    Jets _recojets;
    size_t _jetEffFnHash, _bTagEffFnHash, _cTagEffFnHash, _jetSmearFnHash;
    std::function<double(const Jet&)> _jetEffFn, _bTagEffFn, _cTagEffFn;
    std::function<Jet(const Jet&)> _jetSmearFn;

  };


}

#endif
