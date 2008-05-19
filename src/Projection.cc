// -*- C++ -*-
#include "Rivet/Projection.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Event.hh"
#include "Rivet/Cmp.hh"

namespace Rivet {


  Projection::Projection()
    : _name("BaseProjection")
  {
    addBeamPair(ANY, ANY);
    getLog() << Log::TRACE << "Creating " << getName() << " at " << this << endl;
  }


  Projection:: ~Projection() {
    getLog() << Log::TRACE << "Destroying " << getName() << " at " << this << endl;
  }


  int Projection::compare(const Projection& p) const {
    return mkNamedPCmp(p, "FS");
  }


  bool Projection::before(const Projection& p) const {
    const std::type_info& thisid = typeid(*this);
    const std::type_info& otherid = typeid(p);
    if (thisid == otherid) {
      return compare(p) < 0;
    } else {
      return thisid.before(otherid);
    }
  }
  
  
  const Cuts Projection::getCuts() const {
    Cuts totalCuts = _cuts;
    set<ConstProjectionPtr> projs = getProjections();
    for (set<ConstProjectionPtr>::const_iterator p = projs.begin(); p != projs.end(); ++p) {
      totalCuts.addCuts((*p)->getCuts());
    }
    return totalCuts;
  }
  
  
  const set<BeamPair> Projection::getBeamPairs() const {
    set<BeamPair> ret = _beamPairs;
    set<ConstProjectionPtr> projs = getProjections();
    for (set<ConstProjectionPtr>::const_iterator ip = projs.begin(); ip != projs.end(); ++ip) {
      ConstProjectionPtr p = *ip;
      getLog() << Log::TRACE << "Proj addr = " << p << endl;
      if (p) ret = intersection(ret, p->getBeamPairs());
    }
    return ret;
  }


  Cmp<Projection> Projection::mkNamedPCmp(const Projection& otherparent, 
                                          const string& pname) const {
    return pcmp(*this, otherparent, pname);
  }


  Cmp<Projection> Projection::mkPCmp(const Projection& otherparent, 
                                     const string& pname) const {
    return pcmp(*this, otherparent, pname);
  }


}
