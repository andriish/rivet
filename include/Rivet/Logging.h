// $Id: $
#ifndef RIVET_LOGGING_H 
#define RIVET_LOGGING_H 1

#include "log4cpp/Category.hh"
#include "log4cpp/CategoryStream.hh"
#include "log4cpp/FileAppender.hh"
#include "log4cpp/OstreamAppender.hh"
#include "log4cpp/BasicLayout.hh"
#include <string>

namespace Rivet {
  typedef log4cpp::Category Logger;
  typedef log4cpp::Priority LogPriority;
  static log4cpp::CategoryStream::Separator endlog = log4cpp::CategoryStream::ENDLINE;

  Logger& getLogger(std::string logfile = "");
}

#endif // RIVET_LOGGING_H
