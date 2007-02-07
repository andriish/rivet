// $Id: $
#ifndef RIVET_LOGGING_H 
#define RIVET_LOGGING_H 1

#include "log4cpp/Category.hh"
#include "log4cpp/CategoryStream.hh"
#include "log4cpp/FileAppender.hh"
#include "log4cpp/OstreamAppender.hh"
#include "log4cpp/BasicLayout.hh"
#include <string>

#include "Rivet/Analysis/Analysis.fhh"

namespace Rivet {
  typedef log4cpp::Category Logger;
  typedef log4cpp::Priority LogPriority;
  typedef log4cpp::CategoryStream LogStream;
  static const LogStream::Separator endlog = LogStream::ENDLINE;

  Logger& getLogger(std::string logfile = "");
  Logger& getLogger(AnalysisName analysis);
}

#endif // RIVET_LOGGING_H
