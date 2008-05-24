// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/ProjectionHandler.hh"
#include "Rivet/Cmp.hh"
#include <algorithm>

namespace Rivet {


  // Initialize instance pointer to null.
  ProjectionHandler* ProjectionHandler::_instance = 0;


  ProjectionHandler* ProjectionHandler::create() {
    if (!_instance) {
      _instance = new ProjectionHandler();
      Log::getLog("Rivet.ProjectionHandler") << Log::TRACE << "Created new ProjectionHandler at " 
                                       << _instance << endl;
    }
    return _instance;
  }


  // Get a logger.
  Log& ProjectionHandler::getLog() const {
    return Log::getLog("Rivet.ProjectionHandler");
  }


  void ProjectionHandler::clear() {
    for (ProjHandles::iterator ph = _projs.begin(); ph != _projs.end(); ++ph) {
      getLog() << Log::TRACE << "Deleting projection at " << *ph << endl;
      delete *ph;
    }
    _projs.clear();
    _namedprojs.clear();
  }

  
  // Delete contained pointers.
  ProjectionHandler::~ProjectionHandler() {
    clear();
  }


  void ProjectionHandler::declareClone(const Projection* oldproj, const Projection* newproj) {
    if (oldproj == newproj) return;
    NamedProjsMap::const_iterator nps = _namedprojs.find(oldproj);
    if (nps != _namedprojs.end()) {
      getLog() << Log::TRACE << "Cloning registered projections list: " 
               << oldproj << " -> " << newproj << endl;
      _namedprojs[newproj] = nps->second;
    }
  }


  // Take a {@c clone}'d Projection, compare it to the others on record, and
  // return (by reference) an equivalent Projection which is guaranteed to be
  // the version that will be applied to an event.
  const Projection& ProjectionHandler::registerClonedProjection(const ProjectionApplier& parent, 
                                                                const Projection& oldproj, 
                                                                const Projection& newproj, 
                                                                const string& name) {
    declareClone(&oldproj, &newproj);
    return registerProjection(parent, newproj, name);
  }


  // Take a Projection, compare it to the others on record, and
  // return (by reference) an equivalent Projection which is guaranteed to be
  // the version that will be applied to an event.
  const Projection& ProjectionHandler::registerProjection(const ProjectionApplier& parent, 
                                                          const Projection& proj, 
                                                          const string& name) {
    getLog() << Log::TRACE << "Trying to register projection " << &proj 
             << " for parent " << &parent << " with name '" << name << "'" << endl;

    // Try to find an exact match by pointer
    if (find(_projs.begin(), _projs.end(), &proj) != _projs.end()) {
      getLog() << Log::TRACE << "Matched " << &proj 
               << " by pointer: no cmp loop needed" << endl;
      _namedprojs[&parent][name] = &proj;
      return proj;
    }

    // Get class type using RTTI
    const std::type_info& newtype = typeid(proj);
    getLog() << Log::TRACE << "RTTI type of " << &proj << " is " << newtype.name() << endl; 

    // Compare to ALL projections (use caching _projs).
    getLog() << Log::TRACE << "Comparing " << &proj 
             << " with " << _projs.size()
             << " registered projections" <<  endl;
    for (ProjHandles::const_iterator ph = _projs.begin(); ph != _projs.end(); ++ph) {
      // Make sure the concrete types match, using RTTI.
      const std::type_info& regtype = typeid(**ph);
      getLog() << Log::TRACE << "RTTI type comparison with "<< *ph << ": " 
               << newtype.name() << " vs. " << regtype.name() << endl; 
      if (newtype != regtype) continue;
      getLog() << Log::TRACE << "RTTI type matches with " << *ph << endl;
      
      // If we find a match, ~~delete the passed object, then~~ make a copy of the
      // existing pointer to the store location indexed by ProjApplier* => name.
      if (pcmp(**ph, proj) == PCmp::EQUIVALENT) {
        getLog() << Log::TRACE << "Deleting equivalent projection at " 
                 << &proj << " and returning " << *ph << endl;
        //delete &proj;
        _namedprojs[&parent][name] = *ph;
        return **ph;
      }
    }

    // If we found no match, add passed Projection to _projs and the
    // ProjApplier* => name location in the associative container.
    getLog() << Log::TRACE << "Registered new projection at " << &proj << endl;
    _projs.push_back(&proj);
    _namedprojs[&parent][name] = &proj;
    return proj;
  }


  set<const Projection*> ProjectionHandler::getChildProjections(const ProjectionApplier& parent,
                                                                ProjDepth depth) const {
    set<const Projection*> toplevel;
    NamedProjs nps = _namedprojs.find(&parent)->second;
    for (NamedProjs::const_iterator np = nps.begin(); np != nps.end(); ++np) {
      toplevel.insert(np->second);
    }
    if (depth == SHALLOW) {
      // Only return the projections directly contained within the top level
      return toplevel;
    } else {
      // Return recursively built projection list
      set<const Projection*> alllevels = toplevel;
      for (set<const Projection*>::const_iterator pp = toplevel.begin(); pp != toplevel.end(); ++pp) {
        set<const Projection*> allsublevels = getChildProjections(**pp, DEEP);
        alllevels.insert(allsublevels.begin(), allsublevels.end());
        break;
      }
      return alllevels;
    }
  }



  const Projection& ProjectionHandler::getProjection(const ProjectionApplier& parent,
                                                     const string& name) const {
    getLog() << Log::TRACE << "Searching for child projection '" 
             << name << "' of " << &parent << endl;
    NamedProjsMap::const_iterator nps = _namedprojs.find(&parent);
    if (nps == _namedprojs.end()) {
      ostringstream msg;
      msg << "No projections registered for parent " << &parent;
      throw Error(msg.str());
    }
    NamedProjs::const_iterator np = nps->second.find(name);
    if (np == nps->second.end()) {
      ostringstream msg;
      msg << "No projection '" << name << "' found for parent " << &parent;
      throw Error(msg.str());
    }
    // If it's registered with the projection handler, we must be able to safely
    // dereference the Projection pointer to a reference...
    return *(np->second);
  }



}
