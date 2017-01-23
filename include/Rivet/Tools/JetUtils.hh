#ifndef RIVET_JETUTILS_HH
#define RIVET_JETUTILS_HH

#include "Rivet/Jet.hh"
#include "Rivet/Tools/ParticleBaseUtils.hh"

namespace Rivet {


  /// std::function instantiation for functors taking a Jet and returning a bool
  using JetSelector = function<bool(const Jet&)>;
  /// std::function instantiation for functors taking two Jets and returning a bool
  using JetSorter = function<bool(const Jet&, const Jet&)>;


  /// @name Unbound functions for converting between Jets, Particles and PseudoJets
  //@{

  inline PseudoJets mkPseudoJets(const Particles& ps) {
    PseudoJets rtn; rtn.reserve(ps.size());
    for (const Particle& p : ps)
      rtn.push_back(p);
    return rtn;
  }

  inline PseudoJets mkPseudoJets(const Jets& js) {
    PseudoJets rtn; rtn.reserve(js.size());
    for (const Jet& j : js)
      rtn.push_back(j);
    return rtn;
  }

  inline Jets mkJets(const PseudoJets& pjs) {
    Jets rtn; rtn.reserve(pjs.size());
    for (const PseudoJet& pj : pjs)
      rtn.push_back(pj);
    return rtn;
  }

  //@}


  /// @name Unbound functions for filtering jets
  //@{

  /// Filter a jet collection in-place to the subset that passes the supplied Cut
  Jets& ifilter_select(Jets& jets, const Cut& c);
  /// Alias for ifilter_select
  /// @deprecated Use ifilter_select
  inline Jets& ifilterBy(Jets& jets, const Cut& c) { return ifilter_select(jets, c); }

  /// Filter a jet collection in-place to the subset that passes the supplied Cut
  inline Jets filter_select(const Jets& jets, const Cut& c) {
    Jets rtn = jets;
    return ifilter_select(rtn, c);
  }
  /// Alias for ifilter_select
  /// @deprecated Use filter_select
  inline Jets filterBy(const Jets& jets, const Cut& c) { return filter_select(jets, c); }

  /// Filter a jet collection in-place to the subset that passes the supplied Cut
  inline Jets filter_select(const Jets& jets, const Cut& c, Jets& out) {
    out = filter_select(jets, c);
    return out;
  }
  /// Alias for ifilter_select
  /// @deprecated Use filter_select
  inline Jets filterBy(const Jets& jets, const Cut& c, Jets& out) { return filter_select(jets, c, out); }


  /// Filter a jet collection in-place to the subset that fails the supplied Cut
  Jets& ifilter_discard(Jets& jets, const Cut& c);

  /// Filter a jet collection in-place to the subset that fails the supplied Cut
  inline Jets filter_discard(const Jets& jets, const Cut& c) {
    Jets rtn = jets;
    return ifilter_discard(rtn, c);
  }

  /// Filter a jet collection in-place to the subset that fails the supplied Cut
  inline Jets filter_discard(const Jets& jets, const Cut& c, Jets& out) {
    out = filter_discard(jets, c);
    return out;
  }

  //@}


}

#endif
