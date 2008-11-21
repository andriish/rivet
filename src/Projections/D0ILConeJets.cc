// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Cmp.hh"

namespace Rivet {


  int D0ILConeJets::compare(const Projection& p) const {
    const D0ILConeJets& other = dynamic_cast<const D0ILConeJets&>(p);
    return 
      mkNamedPCmp(other, "FS") || 
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
  }


  void D0ILConeJets::project(const Event& e) {
    // Project into final state
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    calc(fs.particles());
  }


  void D0ILConeJets::reset() { 
    _jets.clear(); 
    _lorentzvecjets.clear();
    _particlelist.clear();
    _particlepointerlist.clear();
  }


  void D0ILConeJets::calc(const ParticleVector& ps) {
    // Store 4 vector data about each particle into list
    foreach (const Particle& p, ps) {
      const HepMC::FourVector fv = p.momentum();
      const HepEntity listelement(fv.e(), fv.px(), fv.py(), fv.pz());
      /// @todo Why both lists?
      _particlelist.push_back(listelement);
      _particlepointerlist.push_back(&_particlelist.back());
    }

    const double ITEM_ET_THRESHOLD = 0.0;
    // Clear jets of previous event
    _jets.clear(); 
    _algo.makeClusters(_jets, _particlepointerlist, ITEM_ET_THRESHOLD); // Compute...
    /// @todo Needed? Does this ordering make sense?
    _particlelist.clear(); // Clear this event
    _particlepointerlist.clear(); // Clear this event
    
    _lorentzvecjets.clear();
    foreach (const HepEntity& jt, _jets) {
      /// @todo Check that ordering is good. Eliminate HepEntity
      FourMomentum jet(jt.E, jt.px, jt.py, jt.pz);
      _lorentzvecjets.push_back(jet);
    }
  }



  Jets D0ILConeJets::jets(double ptmin) const {
    Jets rtn;
    foreach (const HepEntity& he, _jets) {
      // getLog() << Log::TRACE << "D0 HepEntity = (" 
      //          << he.E/GeV << "; " 
      //          << he.px/GeV << ", " 
      //          << he.py/GeV << ", " 
      //          << he.pz/GeV << endl;
      FourMomentum v4(he.E, he.px, he.py, he.pz);
      //getLog() << Log::TRACE << "D0 HepEntity -> V4 = " << v4 << endl; 
      if (v4.pT() < ptmin) continue;
      assert(v4.pT() > ptmin);
      Jet j;
      j.addParticle(v4);
      getLog() << Log::TRACE << "D0 jet pT = " << j.momentum()/GeV << endl;
      rtn.push_back(j);
    }
    // foreach (const Jet& j, rtn) {
    //   getLog() << Log::TRACE << "D0 HepEntity -> V4 -> Jet = " << j.momentum() << endl; 
    // }
    return rtn;
  }
  
  
  void D0ILConeJets::_init_extra_params() {
    _split_ratio = 0.5; 
    _far_ratio = 0.5;
    _et_min_ratio = 0.5; 
    _kill_duplicate = true; 
    _duplicate_dR = 0.005;
    _duplicate_dPT = 0.01; 
    _search_factor = 1.0; 
    _pT_min_leading_protojet = 0.0; 
    _pT_min_second_protojet = 0.0; 
    _merge_max = 1000;
    _pT_min_nomerge = 0.0; 
    _algo = ILConeAlgorithm<HepEntity>(_cone_radius, 
      _min_jet_Et, _split_ratio,_far_ratio, _et_min_ratio, 
      _kill_duplicate, _duplicate_dR, _duplicate_dPT, 
      _search_factor, _pT_min_leading_protojet, 
      _pT_min_second_protojet,_merge_max, _pT_min_nomerge);
  }
  

}
