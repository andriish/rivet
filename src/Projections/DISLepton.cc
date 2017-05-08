// -*- C++ -*-
#include "Rivet/Projections/DISLepton.hh"

namespace Rivet {


  int DISLepton::compare(const Projection& p) const {
    const DISLepton& other = pcast<DISLepton>(p);
    return mkNamedPCmp(other, "Beam") || mkNamedPCmp(other, "FS");
  }


  void DISLepton::project(const Event& e) {
    const ParticlePair& inc = applyProjection<Beam>(e, "Beam").beams();

    bool firstIsLepton = PID::isLepton(inc.first.pid());
    bool secondIsLepton = PID::isLepton(inc.second.pid());

    if (firstIsLepton && !secondIsLepton) {
      _incoming = inc.first;
    } else if (!firstIsLepton && secondIsLepton) {
      _incoming = inc.second;
    } else {
      //eek!
      throw Error("DISLepton projector could not find the correct beam.");
    }
    const GenParticle* current_l=_incoming.genParticle();
    bool found_next_vertex = true;
    while (found_next_vertex) {
        found_next_vertex = false;
        if (!current_l->end_vertex()) break;
        std::vector<const GenParticle*> out_n;
        std::vector<const GenParticle*> out_c;
        for (const GenParticle* pp : particles_out(current_l, HepMC::children)) {
            if (current_l->pdg_id() == pp->pdg_id()) out_n.push_back(pp);
            //+-1 should allow neutrino to electron and electron to neutrino
            if (std::abs(std::abs(current_l->pdg_id()) - std::abs(pp->pdg_id())) == 1) out_c.push_back(pp);
          }
        if (out_n.empty() && out_c.empty()) {
          MSG_WARNING("DISLepton projector: no electron/lepton in the new vertex.");
          break;
        }
        if (out_c.size() + out_n.size() > 1) {
          MSG_WARNING("DISLepton projector: more than one electron/lepton in the new vertex.");
          break;
        }
        if (out_c.size() == 1) current_l = out_c.front();
        if (out_n.size() == 1) current_l = out_n.front();
        found_next_vertex = true;
      }
    if (current_l == NULL)
      throw Error("DISLepton projector could not find the scattered lepton.");
    _outgoing = Particle(current_l);
    if (_outgoing.charge() == _incoming.charge()) _charged=0;
    else _charged = 0; // We consider only electric charge
  }


}
