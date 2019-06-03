// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/ProjectionHandler.hh"
#include "Rivet/Tools/Cmp.hh"
#include <algorithm>
#include <iostream>

using std::cerr;

namespace {
  // Get a logger.
  Rivet::Log& getLog() {
    return Rivet::Log::getLog("Rivet.ProjectionHandler");
  }
}

namespace Rivet {


  // Take a Projection, compare it to the others on record, and return (by
  // reference) an equivalent Projection which is guaranteed to be the
  // (persistent) version that will be applied to an event.
  const Projection& ProjectionHandler::registerProjection(const ProjectionApplier& parent,
                                                          const Projection& proj,
                                                          const string& name)
  {
    getLog() << Log::TRACE << "Trying to register"
             << " projection " << &proj  << " (" << proj.name() << ")"
             << " for parent " << &parent << " (" << parent.name() << ")"
             << " with name '" << name << "'" << endl;

    // Check for duplicate use of "name" on "parent"
    const bool dupOk = _checkDuplicate(parent, proj, name);
    if (!dupOk) {
      cerr << "Duplicate name '" << name << "' in parent '" << parent.name() << "'." << endl;
      exit(1);
    }

    // Choose which version of the projection to register with this parent and name
    ProjHandle ph = _getEquiv(proj);
    if ( ph ) {
      const Projection & ret = _register(parent, ph, name);
      return ret;
    } else {
      unique_ptr<Projection> p = _clone(proj);
      const Projection & ret = _register(parent, move(p), name);
      // Return registered proj
      return ret;
    }
  }

  // Clone neatly
  unique_ptr<Projection> ProjectionHandler::_clone(const Projection& proj)
  {
    // Clone a new copy of the passed projection on the heap
    getLog() << Log::TRACE << "Cloning projection " << proj.name() << " from " << &proj << "..." << endl;
    unique_ptr<Projection> newproj = proj.clone();
    getLog() << Log::TRACE << "...cloned to " << proj.name() << " at " << newproj.get() << endl;

    // Copy all the child ProjHandles when cloning, since otherwise links to "stack parents"
    // will be generated by their children, without any connection to the cloned parent
    if (&proj != newproj.get()) {
      auto nps = _namedprojs.find(&proj);
      if (nps != _namedprojs.end()) {
        getLog() << Log::TRACE << "Cloning registered projections list: "
                 << &proj << " -> " << newproj.get() << endl;
        getLog() << Log::TRACE  << "** creates " << newproj.get() << " -> (map from " << nps->first << ")\n";
        _namedprojs[newproj.get()] = nps->second;
      }
    }

    return newproj;
  }



  // Take a Projection, and register it in the registry.
  const Projection& ProjectionHandler::_register(const ProjectionApplier& parent,
                                                 ProjHandle p,
                                                 const string& name)
  {
    // here we take ownership of the projection
    getLog() << Log::TRACE << "Registering new projection at " << p.get()
      << ". Starting refcount: " << p.use_count() << endl;

    // Add the passed Projection to _projs

    _projs.insert(p);
    getLog() << Log::TRACE
      << "** inserted " << p.get() << " to lookup. Refcount: " << p.use_count() << endl;


    // Add the ProjApplier* => name location to the associative container
    _namedprojs[&parent][name] = p;
    getLog() << Log::TRACE
      << "** created " << &parent << " -> (" << name << ',' <<
                                            p.get() << "). Refcount: " << p.use_count() << endl;

    p->markAsOwned();

    return *p;
  }




  // Try to find a equivalent projection in the system
  ProjHandle ProjectionHandler::_getEquiv(const Projection& proj) const
  {
    // Get class type using RTTI
    const std::type_info& newtype = typeid(proj);
    getLog() << Log::TRACE << "RTTI type of " << &proj << " is " << newtype.name() << endl;

    // Compare to ALL projections via _projs collection
    getLog() << Log::TRACE << "Comparing " << &proj
             << " with " << _projs.size()
             << " registered projection" << (_projs.size() == 1 ? "" : "s") <<  endl;
    for (const ProjHandle& ph : _projs) {
      // Make sure the concrete types match, using RTTI.
      const std::type_info& regtype = typeid(*ph);
      getLog() << Log::TRACE << "  RTTI type comparison with " << ph << ": "
               << newtype.name() << " vs. " << regtype.name() << endl;
      if (newtype != regtype) continue;
      getLog() << Log::TRACE << "  RTTI type matches with " << ph << endl;

      // Test for semantic match
      if (pcmp(*ph, proj) != CmpState::EQ) {
        getLog() << Log::TRACE << "  Projections at "
                 << &proj << " and " << ph << " are not equivalent" << endl;
      } else {
        getLog() << Log::TRACE << "  MATCH! Projections at "
                 << &proj << " and " << ph << " are equivalent" << endl;
        return ph;
      }
    }
    getLog() << Log::TRACE << "  Nothing matches." << endl;
    // If no match, just return a null pointer
    return nullptr;
  }



  string ProjectionHandler::_getStatus() const {
    ostringstream msg;
    msg << "Current projection hierarchy:" << endl;
    for (const NamedProjsMap::value_type& nps : _namedprojs) {
      //const string parentname = nps.first->name();
      msg << nps.first << endl; //"(" << parentname << ")" << endl;
      for (const NamedProjs::value_type& np : nps.second) {
        msg << "  " << np.second << " (" << np.second->name()
            << ", locally called '" << np.first << "')" << endl;
      }
      msg << endl;
    }
    return msg.str();
  }



  // Check that the same parent hasn't already used this name for something else
  bool ProjectionHandler::_checkDuplicate(const ProjectionApplier& parent,
                                          const Projection& proj,
                                          const string& name) const
  {
    auto listedParent = _namedprojs.find(&parent);
    if (listedParent != _namedprojs.end()) {
      const NamedProjs pnps = listedParent->second;
      const NamedProjs::const_iterator ipph = pnps.find(name);
      if (ipph != pnps.end()) {
        const ProjHandle pph = ipph->second;
        getLog() << Log::ERROR << "Projection clash! "
                 << parent.name() << " (" << &parent << ") "
                 << "is trying to overwrite its registered '" << name << "' "
                 << "projection (" << pph << "="
                 << pph->name() << ") with a non-equivalent projection "
                 << "(" << &proj << "=" << proj.name() << ")" << endl;
        getLog() << Log::ERROR << _getStatus();
        return false;
      }
    }
    return true;
  }




  void ProjectionHandler::removeProjectionApplier(ProjectionApplier& parent) {
    auto npi = _namedprojs.find(&parent);
    if (npi != _namedprojs.end()) {
       getLog() << Log::TRACE << "REMOVE Projection at "
                 << &parent << " from map" << endl;
      _namedprojs.erase(npi);
    }
    //
    auto pAsProj = dynamic_cast<Projection*>(&parent);
    if (pAsProj) {
      auto pi = find_if(_projs.begin(), _projs.end(),
                        [pAsProj](ProjHandle h)->bool { return h.get() == pAsProj; } );
      if (pi != _projs.end()) {
          getLog() << Log::TRACE << "REMOVE Projection at "
                   << pAsProj << " from lookup" << endl;
          _projs.erase(pi);

      }
    }
  }


  set<const Projection*> ProjectionHandler::getChildProjections(const ProjectionApplier& parent, ProjDepth depth) const {
    set<const Projection*> toplevel;
    NamedProjs nps = _namedprojs.find(&parent)->second;
    for (NamedProjs::value_type& np : nps) {
      toplevel.insert(np.second.get());
    }
    if (depth == SHALLOW) {
      // Only return the projections directly contained within the top level
      return toplevel;
    } else {
      // Return recursively built projection list
      set<const Projection*> alllevels = toplevel;
      for (const Projection* p : toplevel) {
        set<const Projection*> allsublevels = getChildProjections(*p, DEEP);
        alllevels.insert(allsublevels.begin(), allsublevels.end());
      }
      return alllevels;
    }
  }


  bool ProjectionHandler::hasProjection(const ProjectionApplier& parent, const string& name) const {
    MSG_TRACE("Searching for child projection '" << name << "' of " << &parent);
    NamedProjsMap::const_iterator nps = _namedprojs.find(&parent);
    if (nps == _namedprojs.end()) return false;
    NamedProjs::const_iterator np = nps->second.find(name);
    return !(np == nps->second.end());
  }


  const Projection& ProjectionHandler::getProjection(const ProjectionApplier& parent, const string& name) const {
    MSG_TRACE("Searching for child projection '" << name << "' of " << &parent);
    NamedProjsMap::const_iterator nps = _namedprojs.find(&parent);
    if (nps == _namedprojs.end()) {
      std::ostringstream msg;
      msg << "No projections registered for parent " << &parent;
      throw Error(msg.str());
    }
    NamedProjs::const_iterator np = nps->second.find(name);
    if (np == nps->second.end()) {
      std::ostringstream msg;
      msg << "No projection '" << name << "' found for parent " << &parent;
      throw Error(msg.str());
    }
    MSG_TRACE("Found projection '" << name << "' of " << &parent << " -> " << np->second);
    // If it's registered with the projection handler, we must be able to safely
    // dereference the Projection pointer to a reference...
    return *(np->second);
  }



}
