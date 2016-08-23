#ifndef RIVET_JETUTILS_HH
#define RIVET_JETUTILS_HH

#include "Rivet/Jet.hh"
#include "Rivet/Tools/ParticleBaseUtils.hh"

namespace Rivet {


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

  // /// Filter a jet collection in-place to the subset that passes the supplied function
  // template <typename FN>
  // inline Jets& filterBy(Jets& jets, const FN& f) {
  //   /// @todo Need to invert f logic!
  //   const auto newend = std::remove_if(jets.begin(), jets.end(), f);
  //   jets.erase(newend, jets.end());
  //   return jets;
  // }

  // /// Get a subset of the supplied jets that passes the supplied function
  // template <typename FN>
  // inline Jets filterBy(const Jets& jets, const FN& f) {
  //   Jets rtn = jets;
  //   return ifilterBy(rtn, f);
  // }

  // /// Filter a jet collection to the subset that passes the supplied function, into a new container
  // /// @note New container will be replaced, not appended to
  // template <typename FN>
  // inline Jets& filterBy(Jets& jets, const FN& f, Jets& out) {
  //   out = filterBy(jets, f);
  //   return out;
  // }


  /// Filter a jet collection in-place to the subset that passes the supplied Cut
  Jets& ifilterBy(Jets& jets, const Cut& c);

  /// Get a subset of the supplied jets that passes the supplied Cut
  inline Jets filterBy(const Jets& jets, const Cut& c) {
    Jets rtn = jets;
    return ifilterBy(rtn, c);
  }

  /// Filter a jet collection to the subset that passes the supplied Cut, into a new container
  /// @note New container will be replaced, not appended to
  inline Jets& filterBy(Jets& jets, const Cut& c, Jets& out) {
    out = filterBy(jets, c);
    return out;
  }

  //@}


}

#endif
