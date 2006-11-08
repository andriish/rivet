// $Id: $
#ifndef RIVET_LOGGING_H 
#define RIVET_LOGGING_H 1

/// @Todo Wouldn't it be nice if there was a simple fallback method 
/// that didn't use log4cpp but used the same interface
#include "log4cpp/Category.hh"
#include "log4cpp/FileAppender.hh"
#include "log4cpp/OstreamAppender.hh"
#include "log4cpp/BasicLayout.hh"
#include <string>

namespace Rivet {
  typedef log4cpp::Category Logger;

  Logger& getLogger(std::string logfile = "") {
    static bool ready(false);
    static log4cpp::Layout* layout = 0;
    static log4cpp::Appender* app = 0;

    Logger& log = Logger::getInstance("rivet");
    log.setAdditivity(false);
    if (!layout) { layout = new log4cpp::BasicLayout(); }
    if (app) delete app;

    if (!logfile.empty()) {
      app = new log4cpp::FileAppender("FileAppender", logfile);
    } else {
      app = new log4cpp::OstreamAppender("StreamAppender", &std::cout);
    }
    app->setLayout(layout);
    log.setAppender(app);
    return log;
  }
  
}
#endif // RIVET_LOGGING_H
