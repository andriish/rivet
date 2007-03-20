// $Id: $

#include <string>
#include <iostream>
#include <ctime>
using namespace std;


enum LogLevel { 
  TRACE = 0, DEBUG = 10, INFO = 20, WARN = 30, ERROR = 40 
};


class Logger {

  friend ostream& operator<<(Logger& log, const LogLevel& level);

public:
  Logger(const string& name) 
    : _name(name), _level(INFO), _nostream(new ostream(0)), _writeTime(true) { }

  Logger(const string& name, const LogLevel& level) 
    : _name(name), _level(level), _nostream(new ostream(0)), _writeTime(true) { }
  
public:
  LogLevel getLevel() const {
    return _level;
  }

  Logger& setLevel(const LogLevel& level) {
    _level = level;
    return *this;
  }

  static string getLevelName(const LogLevel& level) {
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

  string getName() const {
    return _name;
  }

  Logger& setName(const string& name) {
    _name = name;
    return *this;
  }

  bool writeTime() const {
    return _writeTime;
  }

  Logger& writeTime(const bool writeTime) {
    _writeTime = writeTime;
    return *this;
  }

  bool isActive(const LogLevel& level) const {
    return (level >= _level);
  }

private:
  string _name;

  LogLevel _level;

  bool _writeTime;

  ostream* const _nostream;

};
  
  
ostream& operator<<(Logger& log, const LogLevel& level) {
  if (log.isActive(level)) {
    time_t rawtime;
    time(&rawtime);
    char* timestr = ctime(&rawtime);
    timestr[24] = ' ';

    cout << Logger::getLevelName(level) << " ";
    cout << log.getName() << ": ";
    if (log.writeTime()) {
      cout << timestr << " ";
    }
    return cout;
  } else {
    return *(log._nostream);
  }
}


int main() {
  Logger log("foo");
  log << INFO << "hello" << endl;
  log << DEBUG << "hi" << endl;
  log << WARN << "hola" << endl;
  log << DEBUG << "hey" << endl;
  return EXIT_SUCCESS;
}
