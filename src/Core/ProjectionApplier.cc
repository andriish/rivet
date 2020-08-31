// -*- C++ -*-
#include "Rivet/ProjectionApplier.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Event.hh"
#include <iostream>

namespace Rivet {


  // NB. Allow proj registration in constructor by default -- explicitly disable for Analysis
  ProjectionApplier::ProjectionApplier()
    : _allowProjReg(true), _owned(false),
      _projhandler(nullptr) // @@@AK
  {  }


  ProjectionApplier::~ProjectionApplier() {
    if ( ! _owned )
      getProjHandler().removeProjectionApplier(*this);
  }


  const Projection& ProjectionApplier::_applyProjection(const Event& evt,
                                                        const string& name) const {
    const Projection& proj = getProjection(name);
    // cout << "Found projection " << &proj << " -> applying" << '\n';
    return _applyProjection(evt, proj);
  }


  const Projection& ProjectionApplier::_applyProjection(const Event& evt,
                                                        const Projection& proj) const {
    return evt.applyProjection(proj);
  }


  const Projection& ProjectionApplier::_declareProjection(const Projection& proj,
                                                          const string& name) {
    if (!_allowProjReg) {
      std::cerr << "Trying to register projection '"
           << proj.name() << "' outside init phase in '" << this->name() << "'.\n";
      exit(2);
    }
    const Projection& reg = getProjHandler().registerProjection(*this, proj, name);
    return reg;
  }

  void ProjectionApplier::setProjectionHandler(ProjectionHandler& projectionHandler) {
      /// Problem with reference reassignment: see comment below on _projhandler member declaration
      //_projhandler = projHandler;
      _projhandler = &projectionHandler;
      /// @todo AB: Move this into a _syncDeclQueue function to be called both by setProjectionHandler() and declare() - [AK] Why we need to call below from declare?
      while (!_declQueue.empty()) {
        /// @todo AB: only use auto when the type is genuinely awkward to express: here it's just - done
        pair<Projection*, string> obj = _declQueue.front();
        /// @todo AB: should the order be switched, to set the PH on the Proj
        /// *about* to be declared first? That way the setting will cascade up
        /// from deepest level to top-level, as currently. Maybe safer that way?
        const Projection& ret = declareProjection(*(obj.first), obj.second);
        //ret.setProjectionHandler(projectionHandler);
        _declQueue.pop();
      }
    }


}
