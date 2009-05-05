#include "Rivet/Event.hh"
#include "Rivet/Tools/Logging.hh"
#include "HepMC/GenEvent.h"

namespace Rivet {


  // Convert the GenEvent to use conventional units (GeV, mm)
  void _geNormUnits(GenEvent& ge) {
    // Specify units if supported
    #ifdef HEPMC_HAS_UNITS
    ge.use_units(HepMC::Units::GEV, HepMC::Units::MM);
    #endif
  }


  void _geRot180x(GenEvent& ge) {
    for (HepMC::GenEvent::particle_iterator ip = ge.particles_begin(); ip != ge.particles_end(); ++ip) {
      const HepMC::FourVector& mom = (*ip)->momentum();
      (*ip)->set_momentum(HepMC::FourVector(mom.px(), -mom.py(), -mom.pz(), mom.e()));
    }
    for (HepMC::GenEvent::vertex_iterator iv = ge.vertices_begin(); iv != ge.vertices_end(); ++iv) {
      const HepMC::FourVector& pos = (*iv)->position();
      (*iv)->set_position(HepMC::FourVector(pos.x(), -pos.y(), -pos.z(), pos.t()));
    }
  }


  // Convert the GenEvent to use conventional alignment 
  // (proton or electron on +ve z-axis?)
  // For example, FHerwig only produces DIS events in the 
  // unconventional orientation and has to be corrected
  void _geNormAlignment(GenEvent& ge) {
    /// @todo Choose when to do in-place rotation of GE
    if (false) {
      Log::getLog("Event") << Log::TRACE << "Rotating event" << endl;
      _geRot180x(ge);
    }
  }


  Event::Event(const GenEvent& ge)
    : _genEvent(ge), _weight(1.0) 
  {
    // Set the weight if there is one, otherwise default to 1.0
    if (!_genEvent.weights().empty()) {
      _weight = ge.weights()[0];
    }

    // Use Rivet's preferred units if possible
    _geNormUnits(_genEvent);
    
    // Use the conventional alignment
    _geNormAlignment(_genEvent);

    // Debug printout to check that copying/magling has worked
    //_genEvent.print();
  }


  Event::Event(const Event& e)
    : _genEvent(e._genEvent), 
      _weight(e._weight) 
  { 
    //
  }


  Event& Event::operator=(const Event& e) {
    _genEvent = e._genEvent;
    _weight = e._weight;
    return *this;
  }


}
