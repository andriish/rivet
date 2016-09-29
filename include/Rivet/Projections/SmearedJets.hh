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
    /// The jet algorithm is mandatory; the smearing and tagging functions are optional
    SmearedJets(const JetAlg& ja,
                const std::vector<std::function<pair<double,Jet>(const pair<double,Jet>&)> > jetDetFn)
      : SmearedJets(ja, jetDetFn, JET_BTAG_PERFECT, JET_CTAG_PERFECT)
    {    }

    /// @brief Constructor with efficiency and smearing function args
    /// The jet algorithm is mandatory; the smearing and tagging functions are optional
    template <typename J2DFN>
    SmearedJets(const JetAlg& ja,
                const std::vector<std::function<pair<double,Jet>(const pair<double,Jet>&)> > jetDetFn,
                const J2DFN& bTagEffFn)
      : SmearedJets(ja, jetDetFn, bTagEffFn, JET_CTAG_PERFECT)
    {    }

    /// @brief Constructor with efficiency and smearing function args
    /// The jet algorithm is mandatory; the smearing and tagging functions are optional
    template <typename J2DFNa, typename J2DFNb>
    SmearedJets(const JetAlg& ja,
                const std::vector<std::function<pair<double,Jet>(const pair<double,Jet>&)> > jetDetFn,
                const J2DFNa& bTagEffFn,
                const J2DFNb& cTagEffFn)
      : _bTagEffFn(bTagEffFn), _cTagEffFn(cTagEffFn), _jetDetFn(jetDetFn)
    {
      setName("SmearedJets");
      addProjection(ja, "TruthJets");
    }

    /// @brief Constructor with efficiency and smearing function args
    /// The jet algorithm is mandatory; the smearing and tagging functions are optional
    template<typename DJ2DJFN>
    SmearedJets(const JetAlg& ja,
                const DJ2DJFN jetDetFn)
      : SmearedJets(ja, jetDetFn, JET_BTAG_PERFECT, JET_CTAG_PERFECT)
    {    }

    /// @brief Constructor with efficiency and smearing function args
    /// The jet algorithm is mandatory; the smearing and tagging functions are optional
    template <typename DJ2DJFN, typename J2DFN>
    SmearedJets(const JetAlg& ja,
                const DJ2DJFN jetDetFn,
                const J2DFN& bTagEffFn)
      : SmearedJets(ja, jetDetFn, bTagEffFn, JET_CTAG_PERFECT)
    {    }


    /// @brief Constructor with efficiency and smearing function args
    /// The jet algorithm is mandatory; the smearing and tagging functions are optional
    template <typename DJ2DJFN, typename J2DFNa, typename J2DFNb>
    SmearedJets(const JetAlg& ja,
                const DJ2DJFN jetDetFn,
                const J2DFNa& bTagEffFn,
                const J2DFNb& cTagEffFn)
      : _bTagEffFn(bTagEffFn), _cTagEffFn(cTagEffFn)
    {
      _jetDetFn.resize(1);
      _jetDetFn[0] = jetDetFn;
      setName("SmearedJets");
      addProjection(ja, "TruthJets");
    }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(SmearedJets);

    //@}


    /// Compare to another SmearedJets
    int compare(const Projection& p) const {
      const SmearedJets& other = dynamic_cast<const SmearedJets&>(p);
      if (get_address(_jetDetFn[0]) == 0) return UNDEFINED;
      if (get_address(_bTagEffFn) == 0) return UNDEFINED;
      if (get_address(_cTagEffFn) == 0) return UNDEFINED;
      MSG_TRACE("Det hashes = ");
      for(size_t i = 0; i < _jetDetFn.size(); ++i) MSG_TRACE( get_address(_jetDetFn[0]) << "," << get_address(other._jetDetFn[0]) << "; ");
      MSG_TRACE("b-tag hashes = " << get_address(_bTagEffFn) << "," << get_address(other._bTagEffFn) << "; " <<
                "c-tag hashes = " << get_address(_cTagEffFn) << "," << get_address(other._cTagEffFn));
      Cmp<unsigned long> ret = mkPCmp(other, "TruthJets") ||
        cmp(get_address(_bTagEffFn), get_address(other._bTagEffFn)) ||
        cmp(get_address(_cTagEffFn), get_address(other._cTagEffFn));
      for(size_t i = 0; i < _jetDetFn.size(); ++i) ret = ret || cmp(get_address(_jetDetFn[i]), get_address(other._jetDetFn[i]));
      return ret;
    }


    /// Perform the jet finding & smearing calculation
    void project(const Event& e) {
      // Copying and filtering
      const Jets& truthjets = apply<JetAlg>(e, "TruthJets").jetsByPt();
      _recojets.clear(); _recojets.reserve(truthjets.size());
      for (const Jet& j : truthjets) {
        // Efficiency sampling
        pair<double, Jet> tmp = make_pair(1., j);
        for(std::function<pair<double,Jet>(const pair<double,Jet>&)> fn : _jetDetFn) tmp = fn(tmp);
        if(tmp.first <= 0) continue;
        if ( !(tmp.first < 1 && rand01() > tmp.first) ) continue;
        // Re-add constituents & tags if (we assume accidentally) they were lost by the smearing function
        if (tmp.second.particles().empty() && !j.particles().empty()) tmp.second.particles() = j.particles();
        if (tmp.second.tags().empty() && !j.tags().empty()) tmp.second.tags() = j.tags();
        _recojets.push_back(tmp.second);
      }

      // Apply tagging efficiencies, using smeared kinematics as input to the tag eff functions
      for (Jet& j : _recojets) {
        const double beff = _bTagEffFn ? _bTagEffFn(j) : 1;
        const bool btag = beff == 1 || (beff != 0 && beff < rand01());
        // Remove b-tags if needed, and add a dummy one if needed
        if (!btag && j.bTagged()) j.tags().erase(std::remove_if(j.tags().begin(), j.tags().end(), hasBottom), j.tags().end());
        if (btag && !j.bTagged()) j.tags().push_back(Particle(PID::BQUARK, j.mom()));
        const double ceff = _cTagEffFn ? _cTagEffFn(j) : 1;
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

    /// Stored efficiency functions
    std::function<double(const Jet&)> _bTagEffFn, _cTagEffFn;
    std::vector<std::function<pair<double,Jet>(const pair<double,Jet>&)> > _jetDetFn;

  };


}

#endif
