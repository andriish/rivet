// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Cmp.hh"

namespace Rivet {

  /// @todo The compare method is incomplete?
  int D0ILConeJets::compare(const Projection& p) const {
    const D0ILConeJets& other = dynamic_cast<const D0ILConeJets&>(p);
    return 
      pcmp(_fsproj, other._fsproj) || 
      cmp(_cone_radius, other._cone_radius) ||
      cmp(_min_jet_Et, other._min_jet_Et) || 
      cmp(_split_ratio, other._split_ratio) ||   
      cmp(_far_ratio, other._far_ratio) ||
      cmp(_et_min_ratio, other._et_min_ratio) ||
      cmp(_kill_duplicate, other._kill_duplicate) ||
      cmp(_duplicate_dR, other._duplicate_dR) || 
      cmp(_duplicate_dPT, other._duplicate_dPT) || 
      cmp(_search_factor, other._search_factor) || 
      cmp(_pT_min_leading_protojet, other._pT_min_leading_protojet) || 
      cmp(_pT_min_second_protojet, other._pT_min_second_protojet) ||
      cmp(_merge_max, other._merge_max) || 
      cmp(_pT_min_nomerge, other._pT_min_nomerge);
      // || cmp(_jets, other._jets);
  }


  void D0ILConeJets::project(const Event& e) {

    // Project into final state
    const FinalState& fs = e.applyProjection(_fsproj);

    // Store 4 vector data about each particle into list
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      HepMC::FourVector fv = p->getMomentum();
      const HepEntity listelement(fv.e(), fv.px(), fv.py(), fv.pz());
      _particlelist.push_back(listelement);
      _particlepointerlist.push_back(&_particlelist.back());
    }



    float item_ET_Threshold = 0.0;
    clearJets(); //Clear jets of previous event
    _algo.makeClusters(getJets(), _particlepointerlist, item_ET_Threshold); // Turn the crank!!!
    _particlelist.clear(); //Clear this event
    _particlepointerlist.clear(); //Clear this event
    
    _lorentzvecjets.clear();
    for (list<HepEntity>::const_iterator jt = _jets.begin(); jt != _jets.end(); ++jt) {
      /// @todo Check that ordering is good. Eliminate HepEntity
      FourMomentum jet(jt->E, jt->px, jt->py, jt->pz);
      _lorentzvecjets.push_back(jet);
    }

  }


}
