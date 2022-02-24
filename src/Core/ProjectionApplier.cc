// -*- C++ -*-
#include "Rivet/ProjectionApplier.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Event.hh"
#include <iostream>

namespace Rivet {


  // NB. Allow proj registration in constructor by default -- explicitly disable for Analysis
  // TODO @TP : I have a feeling this behaviour might be slightly deprecated after the thread-safety changes,
  // but I'm not fully sure.
  ProjectionApplier::ProjectionApplier()
    : _allowProjReg(true), _owned(false),
    _projhandler(nullptr)
  {  }


  ProjectionApplier::~ProjectionApplier() {
    //todo @TP: Is owned still needed now that _projhandler is a ptr?
    if ( ! _owned && _projhandler != nullptr){
      getProjHandler().removeProjectionApplier(*this);
    }
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
                                                          const string& name) const {
    MSG_TRACE("Declaring Projection "<<&proj<<" ("<<proj.name()<<") in parent "<<this<<" ("<<this->name()<<")");
    if (!_allowProjReg) {
      std::cerr << "Trying to register projection '"
           << proj.name() << "' outside init phase in '" << this->name() << "'.\n";
      exit(2);
    }
    const Projection& reg = getProjHandler().registerProjection(*this, proj, name);
    _syncDeclQueue();
    return reg;
  }

  void ProjectionApplier::setProjectionHandler(ProjectionHandler& projectionHandler) const {
    _projhandler = &projectionHandler;
    //TODO @TP: I don't think this call is needed anymore?
    _syncDeclQueue();
  }

  void ProjectionApplier::_syncDeclQueue() const {
    MSG_TRACE("Flushing declQueue of ProjectionApplier " << this << " (" << this->name() <<")" << std::endl);
    while (!_declQueue.empty()) {
      pair<std::shared_ptr<Projection>, string> obj = _declQueue.front();
      obj.first->setProjectionHandler(getProjHandler());
      _declQueue.pop_front();
      const Projection& ret = declareProjection(*(obj.first), obj.second);
    }
  }
  
}
