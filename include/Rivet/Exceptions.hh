#ifndef RIVET_EXCEPTIONS_HH 
#define RIVET_EXCEPTIONS_HH

#include <string>
#include <exception>
#include <stdexcept>

namespace Rivet {

  /// Basic unspecialised Rivet exception.
  class Exception : public std::exception { 
  public:
    Exception() : std::exception() {}
  };

  /// @brief Generic runtime Rivet error.
  class Error : public std::runtime_error, public Exception {
  public:
    Error(const std::string& what) : std::runtime_error(what) {} 
  };


  /// Error for e.g. use of invalid bin ranges.
  class RangeError : public std::out_of_range, public Exception {
  public:
    RangeError(const std::string& what) : std::out_of_range(what) {} 
  };


  /// @todo Clarify where this might arise!
  class LogicError : public std::logic_error, public Exception {
  public:
    LogicError(const std::string& what) : std::logic_error(what) {} 
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
