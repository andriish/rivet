// -*- C++ -*-
#ifndef RIVET_RivetInfo_H
#define RIVET_RivetInfo_H
//
// This is the declaration of the RivetInfo class.
//

#include <map>
//#include <multimap>
#include <string>
#include <iostream>
#include "Rivet/Config/Rivet.h"

namespace Rivet {

/**
 * RivetInfo contains information which can be passed from the
 * different Projection and AnalysisBase objects in Rivet to the
 * outside world. The info is in the form of integers or floating
 * point numbers mapped to names given as simple strings. The typical
 * information is about cuts assumed in an analysis which can be used
 * by the event generators to speed up the generation.  The names of
 * the parameters are completely implementation-defined, but there are
 * a number of pre-defined names of specific types. Each parameter can
 * be assigned a type and a description. The type is simply given by a
 * string and is also implementation-defined athough, again, there are
 * some pre-defined types.
 * 
 * The following types of parameters are defined:
 * <ul>
 *
 * <li> \c min means the parameter specifies a minimum value of a
 * parameter. If different objects in the same run specify different
 * values of the same parameter, they will be reduced to just one
 * taking the minimum of all values.
 *
 * <li> \c unique means that if different objects in the same run
 * specify this parameter the values must be the same.
 *
 * <li> \c none is the default type. If different objects in the same
 * run specify different values of the same parameter, they will all
 * be reported.
 *
 * <li> \c limited:par A parameter similar to \c none, the variation
 * in which is limited by another parameter named \c par. If different
 * objects in the same run specify different values, the minimum and
 * maximum may not differ by more than the amount specified by \c
 * par. If \c par is missing or zero this is the same as \c unique.
 *
 * </ul>
 *
 * The pre-defined parameter names are:<ul>
 *
 * <li> \c BeamA and \c BeamB (type: \c unique) the particle id
 * numbers of the incoming beams (A is assumed to be along the
 * positive z-axis).
 *
 * <li> \c EBeamA and \c EBeamB (\c limited:ETolBeamA and \c
 * limited:ETolBeamB) the energy of the incoming beams.
 *
 * <li> \c ETolBeamA and \c ETolBeamB (\c min) The tolerance allowed
 * in the energy of the beams.
 *
 * <li> \c MinQ2 (\c min) The minimum \f$ Q^2 \f$ assumed in a DIS
 * event.
 *
 * </ul>
 */
class RivetInfo {

public:

  /**
   * Typedef for maps of parameter types keyed by the parameter name.
   */
  typedef std::map<std::string,std::string> InfoTMap;

  /**
   * Typedef for maps of parameter descriptions keyed by the
   * parameter name.
   */
  typedef std::multimap<std::string,std::string> InfoMap;

  /**
   * Typedef for maps of integer parameters keyed by the parameter
   * name.
   */
  typedef std::multimap<std::string,long> InfoIMap;

  /**
   * Typedef for maps of floating point parameters keyed by the
   * parameter name.
   */
  typedef std::multimap<std::string,double> InfoDMap;


public:

  /** @name Standard constructors, destructorsand assignment. */
  //@{
  /**
   * The default constructor.
   */
  RivetInfo();

  /**
   * The copy constructor.
   */
  inline RivetInfo(const RivetInfo &);

  /**
   * The destructor.
   */
  virtual ~RivetInfo();

  /**
   * The assignment operator.
   */
  RivetInfo & operator=(const RivetInfo &);

  //@}

public:

  /**
   * Declare an integer parameter to be accessible to the outside
   * world. The \a name and the \a value must be given. If this is not
   * a standard parameter, the \a type and \a description should be
   * given as well.
   */
  void declareParameter(std::string name, long value,
			std::string type = "", std::string description = "");

  /**
   * Declare an floating point parameter to be accessible to the
   * outside world. The \a name and the \a value must be given. If
   * this is not a standard parameter, the \a type and \a description
   * should be given as well.
   */
  void declareParameter(std::string name, double value,
			std::string type = "", std::string description = "");

  /**
   * Add a \a description of the parameter with the given \a name.
   */
  void describeParameter(std::string name, std::string description);

  /**
   * Set the \a type of the parameter with the given \a name.
   */
  void parameterType(std::string name, std::string type);

  /**
   * Append all the parameters from \a inf.
   */
  void append(const RivetInfo & inf);

  /**
   * Append all the parameters from \a inf.
   */
  inline RivetInfo & operator+=(const RivetInfo & inf);

  /**
   * Return a new RivetInfo object with all information in this and
   * \a inf added together.
   */
  inline RivetInfo operator+(const RivetInfo & inf) const;

  /**
   * Check consistency and purge multiple parameter
   * definitions. @throw runtime_error if there were conflicting
   * parameters.
   */
  void check();

  /**
   * Return the value of the parameter with the given \a name. If
   * several values exist, the average value is returned. Returns zero
   * if no parameter with corresponding \a name is found.
   */
  double getFloatParameter(std::string name) const;

  /**
   * Return the value of the parameter with the given \a name. If
   * several values exist, the average value is returned. Returns zero
   * if no parameter with corresponding \a name is found.
   */
  long getIntParameter(std::string name) const;

  /**
   * Return the description of the parameter with the given \a name.
   */
  std::string getDescription(std::string name) const;

  /**
   * Return the type of the parameter with the given \a name.
   */
  inline std::string getType(std::string name) const;

  /**
   * Print the parameters to the given \a stream.
   */
  std::ostream & print(std::ostream & stream) const;

private:

  /**
   * Map the parameter names of this projection to strings
   * representing their types.
   */
  InfoTMap theParTypes;

  /**
   * Map the parameter names of this projection to strings
   * describing them.
   */
  InfoMap theParDescriptions;

  /**
   * Map the names of the integer parameters of this projection to
   * their value.
   */
  InfoIMap theIntPars;

  /**
   * Map the names of the floating point parameters of this projection
   * to their value.
   */
  InfoDMap theFloatPars;

};

inline std::ostream & operator<<(std::ostream & os, const RivetInfo & i) {
  return i.print(os);
}

}

#include "RivetInfo.icc"

#endif /* RIVET_RivetInfo_H */
