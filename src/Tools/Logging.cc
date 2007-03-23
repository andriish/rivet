// $Id: $

#include "Rivet/Tools/Logging.hh"
#include <ctime>

using namespace std;

// Logger& Rivet::getLogger(std::string logfile) {
//   static log4cpp::Layout* layout = 0;
//   static log4cpp::Appender* app = 0;
  
//   if (Logger* p_log = Logger::exists("rivet")) {
//     return *p_log;
//   } else {
//     Logger& log = Logger::getInstance("rivet");
//     log.setAdditivity(false);
//     if (!layout) { layout = new log4cpp::BasicLayout(); }
//     if (app) delete app;
    
//     if (!logfile.empty()) {
//       app = new log4cpp::FileAppender("FileAppender", logfile);
//     } else {
//       app = new log4cpp::OstreamAppender("StreamAppender", &std::cout);
//     }
//     app->setLayout(layout);
//     log.setAppender(app);
//     return log;
//   }
// } 

// /// @todo Get logger priorities by analysis from RivetHandler
// Logger& Rivet::getLogger(AnalysisName analysis) {
//   Logger& log = getLogger();
//   log.setPriority(LogPriority::INFO);
//   return log;
// }

namespace Rivet {

  Log::Log(const string& name) 
    : _name(name), _level(INFO), _writeTime(true), _nostream(new ostream(0)) { }

  Log::Log(const string& name, const Level& level) 
    : _name(name), _level(level), _writeTime(true), _nostream(new ostream(0)) { }

  Log& Log::getLog(const string& name) {
    return *(new Log(name));
  }

  Log& Log::getLog(const string& name, const Level& level) {
    return *(new Log(name, level));
  }

  string Log::getLevelName(const Level& level) {
    switch(level) {
    case TRACE:
      return "TRACE";
    case DEBUG:
      return "DEBUG";
    case INFO:
      return "INFO";
    case WARN:
      return "WARN";
    case ERROR:
      return "ERROR";
    }
  }


  string Log::formatMessage(const Level& level, const std::string& message) {
    string out;
    out += Log::getLevelName(level);
    out += " ";

    out += getName();
    out += ": ";

    if (writeTime()) {
      time_t rawtime;
      time(&rawtime);
      char* timestr = ctime(&rawtime);
      timestr[24] = ' ';
      out += timestr;
      out += " ";
    }

    out += message;
    return out;
  }


  void Log::log(const Level& level, const std::string& message) {
    if (isActive(level)) {
      cout << formatMessage(level, message) << endl;
    }
  }


  ostream& operator<<(Log& log, const Log::Level& level) {
    if (log.isActive(level)) {
      cout << log.formatMessage(level, "");
      return cout;
    } else {
      return *(log._nostream);
    }
  }

}
