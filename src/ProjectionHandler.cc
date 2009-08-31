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
      Log::getLog("Rivet.ProjectionHandler") 
        << Log::TRACE << "Created new ProjectionHandler at " << _instance << endl;
    }
    return _instance;
  }


  // Get a logger.
  Log& ProjectionHandler::getLog() const {
    return Log::getLog("Rivet.ProjectionHandler");
  }


  void ProjectionHandler::clear() {
    foreach (ProjHandles::value_type& ph, _projs) {
      getLog() << Log::TRACE << "Deleting projection at " << ph << endl;
      delete ph;
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
    getLog() << Log::TRACE << "Trying to register"
             << " projection " << &proj  << "(" << proj.name() << ")"
             << " for parent " << &parent << "(" << parent.name() << ")"
             << " with name '" << name << "'" << endl;

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
    foreach (const ProjHandle& ph, _projs) {
      // Make sure the concrete types match, using RTTI.
      const std::type_info& regtype = typeid(*ph);
      getLog() << Log::TRACE << "RTTI type comparison with "<< ph << ": " 
               << newtype.name() << " vs. " << regtype.name() << endl; 
      if (newtype != regtype) continue;
      getLog() << Log::TRACE << "RTTI type matches with " << ph << endl;
      
      // If we find a match, ~~delete the passed object, then~~ make a copy of the
      // existing pointer to the store location indexed by ProjApplier* => name.
      if (pcmp(*ph, proj) != PCmp::EQUIVALENT) {
        getLog() << Log::TRACE << "Type-matched projections at " 
                 << &proj << " and " << ph << " are not equivalent" << endl;
      } else {
        getLog() << Log::TRACE << "Deleting equivalent projection at " 
                 << &proj << " and returning " << ph << endl;
        //delete &proj;
        _namedprojs[&parent][name] = ph;
        return *ph;
      }
    }


    /// @todo Reinstate this check, since something is going wrong with undead projections!
    // // If there is no match, check that the same parent hasn't already used this name for something else
    // if (_namedprojs[&parent].find(name) != _namedprojs[&parent].end()) {
    //   getLog() << Log::ERROR << "Projection clash! "
    //            << parent.name() << " (" << &parent << ") "
    //            << "is trying to overwrite its registered '" << name << "' " 
    //            << "projection (" << _namedprojs[&parent][name] << "=" 
    //            << _namedprojs[&parent][name]->name() << ") with a non-equivalent projection "
    //            << "(" << &proj << "=" << proj.name() << ")" << endl;
    //   ostringstream msg;
    //   msg << "Current projection hierarchy:" << endl;
    //   foreach (const NamedProjsMap::value_type& nps, _namedprojs) {
    //     //const string parentname = nps.first->name();
    //     msg << nps.first << endl; //"(" << parentname << ")" << endl;
    //     foreach (const NamedProjs::value_type& np, nps.second) {
    //       msg << "  " << np.second << " (" << np.second->name() 
    //            << ", locally called '" << np.first << "')" << endl;
    //     }
    //     msg << endl;
    //   }
    //   getLog() << Log::ERROR << msg.str();
    //   exit(1);
    // }


    // If we found no match, and the name is free, add the passed Projection to _projs, and
    // add the ProjApplier* => name location to the associative container.
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
