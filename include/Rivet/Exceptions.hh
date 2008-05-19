#ifndef RIVET_EXCEPTIONS_HH 
#define RIVET_EXCEPTIONS_HH

#include <string>
#include <exception>
#include <stdexcept>

namespace Rivet {

  /// Generic runtime Rivet error.
  class Error : public std::runtime_error {
  public:
    Error(const std::string& what) : std::runtime_error(what) {} 
  };


  /// Also typedef Exception, so that it's there.
  typedef Error Exception;


  /// Error for e.g. use of invalid bin ranges.
  class RangeError : public Error {
  public:
    RangeError(const std::string& what) : Error(what) {} 
  };


  /// @todo Clarify where this might arise!
  class LogicError : public Error {
  public:
    LogicError(const std::string& what) : Error(what) {} 
  };

  /// @brief Errors relating to event/bin weights
  /// Arises in computing statistical quantities because e.g. the bin
  /// weight is zero or negative.
  class WeightError : public Error {
  public:
    WeightError(const std::string& what) : Error(what) {} 
  };

}

#endif
