// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RivetInfo class.
//

#include "RivetInfo.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Rivet;

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

