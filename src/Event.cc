#include "Rivet/Event.hh"
#include "HepMC/GenEvent.h"

namespace Rivet {


  Event::Event(const GenEvent& ge)
    : _genEvent(&ge), _weight(1.0) 
  {
    /// @todo Deep copy of GenEvent (store full GE, not pointer)
    /// @todo Rotation of copy (via external code... don't put the rotation logic here)
    /// @todo Specify units if supported
    if (!ge.weights().empty()) {
      _weight = ge.weights()[0];
    }
  }


  Event::Event(const Event& e)
    : _genEvent(e._genEvent), _weight(e._weight) 
  { }


}
