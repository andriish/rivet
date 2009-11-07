// -*- C++ -*-
#include "Rivet/ProjectionApplier.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Event.hh"

namespace Rivet {


  ProjectionApplier::ProjectionApplier()
    : _projhandler(ProjectionHandler::create())
  { }


  ProjectionApplier::~ProjectionApplier() {
    getProjHandler().removeProjectionApplier(*this);
  }


  const Projection& ProjectionApplier::_applyProjection(const Event& evt, 
                                                        const string& name) const {
    return evt.applyProjection(getProjection(name));
  }
  

  const Projection& ProjectionApplier::_applyProjection(const Event& evt, 
                                                        const Projection& proj) const {
    return evt.applyProjection(proj);
  }


  const Projection& ProjectionApplier::_addProjection(const Projection& proj, 
                                                      const std::string& name) {
    const Projection& reg = getProjHandler().registerProjection(*this, proj, name);
    return reg;
  }


}
