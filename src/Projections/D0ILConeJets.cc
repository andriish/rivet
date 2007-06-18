// -*- C++ -*-

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/Cmp.hh"
#include "Rivet/RivetCLHEP.hh"

#include "Rivet/Tools/D0RunIIcone/inline_maths.h"
#include "Rivet/Tools/D0RunIIcone/HepEntity.h"

namespace Rivet {


  /// @todo The compare method is incomplete?
  int D0ILConeJets::compare(const Projection& p) const {
    const D0ILConeJets& other = dynamic_cast<const D0ILConeJets&>(p);
    return 
      pcmp(*_fsproj, *other._fsproj) || 
      cmp(_cone_radius, other._cone_radius) ||
      cmp(_min_jet_Et, other._min_jet_Et) || 
      cmp(_split_ratio, other._split_ratio);
    //cmp(jets, other.jets);
  }


  void D0ILConeJets::project(const Event& e) {
    /// @todo Enormous memory leak!
    //_particlelist = new list<const HepEntity*>;

    // Project into final state
    const FinalState& fs = e.applyProjection(*_fsproj);

    // Store 4 vector data about each particle into list
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      HepMC::FourVector fv = p->getMomentum();
      /// @todo Why HepEntity* rather than just HepEntity?
      /// @todo Is this ever deleted?
      const HepEntity* listelement = new HepEntity(fv.e(), fv.px(), fv.py(), fv.pz());
      _particlelist->push_back(listelement);
    }

    float item_ET_Threshold = 0.0;
    // jets = pointer to list of type HepEntity
    jets->clear(); // Be sure to have no jets from old event(s)
    _algo->makeClusters(*jets, *_particlelist, item_ET_Threshold); // Turn the crank!!!
    _particlelist->clear();
    
  }


}
