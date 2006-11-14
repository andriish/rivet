// $Id: $

#include "Rivet/Tools/Logging.h"

Rivet::Logger& getLogger(std::string logfile) {
  static log4cpp::Layout* layout = 0;
  static log4cpp::Appender* app = 0;
  
  if (Logger* p_log = Logger::exists("rivet")) {
    return *p_log;
  } else {
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
