// $Id: $

#include "Rivet/Tools/Logging.hh"
#include <ctime>

using namespace std;


namespace Rivet {

  Log::LogMap Log::existingLogs;
  Log::LevelMap Log::defaultLevels;


  Log::Log(const string& name) 
    : _name(name), _level(INFO), _writeTime(true), _nostream(new ostream(0)) { }


  Log::Log(const string& name, const Level& level)
    : _name(name), _level(level), _writeTime(true), _nostream(new ostream(0)) { }


  Log& Log::getLog(const string& name) {
    if (existingLogs.find(name) == existingLogs.end()) {
      Level level = INFO;
      if (defaultLevels.find(name) != defaultLevels.end()) {
        level = defaultLevels.find(name)->second;
      }
      existingLogs[name] = new Log(name, level);
    }
    return *existingLogs[name];
  }


  Log& Log::getLog(const string& name, const Level& level) {
    if (existingLogs.find(name) == existingLogs.end()) {
      existingLogs[name] = new Log(name, level);
    }
    return *existingLogs[name];
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
    throw runtime_error("Enum value was not a valid log level. How did that happen?");
  }


  Log::Level Log::getLevelFromName(const string& level) {
    if (level == "TRACE") return TRACE;
    if (level == "DEBUG") return DEBUG;
    if (level == "INFO") return INFO;
    if (level == "WARN") return WARN;
    if (level == "ERROR") return ERROR;
    throw runtime_error("Couldn't create a log level from string '" + level + "'");
  }


  string Log::formatMessage(const Level& level, const std::string& message) {
    string out;
    out += Log::getLevelName(level);
    out += " ";

    out += getName();
    out += ": ";

    if (isTimeInOutput()) {
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
