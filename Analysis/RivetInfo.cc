// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RivetInfo class.
//

#include "RivetInfo.h"

using namespace Rivet;

RivetInfo::RivetInfo() {
  theParTypes["BeamA"] = "unique";
  theParDescriptions.insert
    (make_pair("BeamA",
	       "The type of the incoming particle along the positive z-axis."));
  theParTypes["BeamB"] = "unique";
  theParDescriptions.insert
    (make_pair("BeamB",
	       "The type of the incoming particle along the negative z-axis."));
}

RivetInfo::~RivetInfo() {}

RivetInfo & RivetInfo::operator=(const RivetInfo & i) {
  theParTypes = i.theParTypes;
  theParDescriptions = i.theParDescriptions;
  theIntPars = i.theIntPars;
  theFloatPars = i.theFloatPars;
  return *this;
}

void RivetInfo::
declareParameter(string name, long value, string type, string description) {
  theIntPars.insert(make_pair(name, value));
  if ( type.length() ) {
    if ( theParTypes.find(name) != theParTypes.end() &&
	 theParTypes[name] != type )
      throw runtime_error("Conflicting types for parameter " + name  +
			  " in RivetInfo object.");
    theParTypes[name] = type;
  }
  if ( description.length() )
    theParDescriptions.insert(make_pair(name, description));
}

void RivetInfo::
declareParameter(string name, double value, string type, string description) {
  theFloatPars.insert(make_pair(name, value));
  if ( type.length() ) {
    if ( theParTypes.find(name) != theParTypes.end() &&
	 theParTypes[name] != type )
      throw runtime_error("Conflicting types for parameter " + name  +
			  " in RivetInfo object.");
    theParTypes[name] = type;
  }
  if ( description.length() )
    theParDescriptions.insert(make_pair(name, description));
}

void RivetInfo::check() {
  set<string> done;
  long oldsize = -1;

  while ( theParTypes.size() < done.size() ) {
    if ( long(done.size()) == oldsize )
      throw runtime_error("Failed to check parameters in RivetInfo object.");
    oldsize = done.size();
    for ( InfoTMap::iterator it = theParTypes.begin();
	  it != theParTypes.end(); ++it ) {
      if ( done.find(it->first) != done.end() ) continue;
      if ( it->second == "unique" ) {
	pair<InfoDMap::iterator,InfoDMap::iterator>
	  range = theFloatPars.equal_range(it->first);
	if ( range.first != theFloatPars.end() ) {
	  double d1 = range.first->second;
	  while ( ++range.first != range.second )
	    if ( d1 != range.first->second )
	      throw runtime_error
		("Found two different values of unique parameter '" +
		 it->first + "'.");
	  theFloatPars.erase(range.first, range.second);
	  theFloatPars.insert(make_pair(it->first, d1));
	} else {
	  pair<InfoIMap::iterator,InfoIMap::iterator>
	    range = theIntPars.equal_range(it->first);
	  if ( range.first != theIntPars.end() ) {
	    long l1 = range.first->second;
	    while ( ++range.first != range.second )
	      if ( l1 != range.first->second )
		throw runtime_error
		  ("Found two different values of unique parameter '" +
		   it->first + "'.");
	    theIntPars.erase(range.first, range.second);
	    theIntPars.insert(make_pair(it->first, l1));
	  }
	}
      }
      else if ( it->second == "min" ) {
	pair<InfoDMap::iterator,InfoDMap::iterator>
	  range = theFloatPars.equal_range(it->first);
	if ( range.first != theFloatPars.end() ) {
	  double dmin = range.first->second;
	  while ( ++range.first != range.second )
	    dmin = min(dmin, range.first->second);
	  theFloatPars.erase(range.first, range.second);
	  theFloatPars.insert(make_pair(it->first, dmin));
	} else {
	  pair<InfoIMap::iterator,InfoIMap::iterator>
	    range = theIntPars.equal_range(it->first);
	  if ( range.first != theIntPars.end() ) {
	    long lmin = range.first->second;
	    while ( ++range.first != range.second )
	      lmin = min(lmin, range.first->second);
	    theIntPars.erase(range.first, range.second);
	    theIntPars.insert(make_pair(it->first, lmin));
	  }
	}
      }
      else if (it->second.substr(0, 8) == "limited:" ) {
	string ref = it->second.substr(8);
	if ( done.find(ref) == done.end() ) {
	  if ( theParTypes.find(ref) == theParTypes.end() )
	    throw runtime_error
	      ("The parameter '" + it->first +
	       "' is limited by the non-existing parameter '" + ref + "'.");
	  continue;
	}
	pair<InfoDMap::iterator,InfoDMap::iterator>
	  range = theFloatPars.equal_range(it->first);
	if ( range.first != theFloatPars.end() ) {
	  double dmin = range.first->second;
	  double dmax = range.first->second;
	  while ( ++range.first != range.second ) {
	    dmin = min(dmin, range.first->second);
	    dmax = min(dmax, range.first->second);
	  }
	  if ( dmin - dmax > getFloatParameter(ref) )
	    throw runtime_error
	      ("The different values of parameter '" + it->first +
	       "' differs by more that is allowed by parameter '" +
	       ref + "'.");
	} else {
	  pair<InfoIMap::iterator,InfoIMap::iterator>
	    range = theIntPars.equal_range(it->first);
	  if ( range.first != theIntPars.end() ) {
	    long lmin = range.first->second;
	    long lmax = range.first->second;
	    while ( ++range.first != range.second ) {
	      lmin = min(lmin, range.first->second);
	      lmax = max(lmax, range.first->second);
	    }
	    if ( lmin - lmax > getIntParameter(ref) )
	      throw runtime_error
		("The different values of parameter '" + it->first +
		 "' differs by more that is allowed by parameter '" +
		 ref + "'.");
	  }
	}
      }
      done.insert(it->first);
    }
  }
}

double RivetInfo::getFloatParameter(string name) const {
  pair<InfoDMap::const_iterator, InfoDMap::const_iterator> range =
    theFloatPars.equal_range(name);
  if ( range.first == theFloatPars.end() ) return 0.0;
  double sum = 0.0;
  double N = 0;
  while ( range.first != range.second ) {
    N += 1.0;
    sum += range.first->second;
    range.first++;
  }
  return sum/N;
}

long RivetInfo::getIntParameter(string name) const {
  pair<InfoIMap::const_iterator, InfoIMap::const_iterator> range =
    theIntPars.equal_range(name);
  if ( range.first == theIntPars.end() ) return 0;
  double sum = 0.0;
  double N = 0;
  while ( range.first != range.second ) {
    N += 1.0;
    sum += range.first->second;
    range.first++;
  }
  return long(sum/N + 0.5);
}


void RivetInfo::append(const RivetInfo & inf) {
  for ( InfoDMap::const_iterator it = inf.theFloatPars.begin();
	it != inf.theFloatPars.end(); ++it )
    declareParameter(it->first, it->second, getType(it->first),
		     getDescription(it->first));
  for ( InfoIMap::const_iterator it = inf.theIntPars.begin();
	it != inf.theIntPars.end(); ++it )
    declareParameter(it->first, it->second, getType(it->first),
		     getDescription(it->first));
  check();
}

