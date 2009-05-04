#include "Rivet/Event.hh"
#include "HepMC/GenEvent.h"

namespace Rivet {


  // Convert the GenEvent to use conventional units (GeV, mm)
  void _geNormUnits(GenEvent& ge) {
    // Specify units if supported
    #ifdef HEPMC_HAS_UNITS
    ge.use_units(HepMC::Units::GEV, HepMC::Units::MM);
    #endif
  }


  // Convert the GenEvent to use conventional alignment 
  // (proton or electron on +ve z-axis?)
  // For example, FHerwig only produces DIS events in the 
  // unconventional orientation and has to be corrected
  void _geNormAlignment(GenEvent& ge) {
    /// @todo In-place rotation of GE
  }


  Event::Event(const GenEvent& ge)
    : _genEvent(ge), _weight(1.0) 
  {
    // Set the weight if there is one, otherwise default to 1.0
    if (!ge.weights().empty()) {
      _weight = ge.weights()[0];
    }

    // Use Rivet's preferred units if possible
    _geNormUnits(_genEvent);
    
    // Use the conventional alignment
    _geNormAlignment(_genEvent);
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
  }


}
