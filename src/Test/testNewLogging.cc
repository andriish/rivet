// $Id: $

#include <string>
#include <iostream>
#include <ctime>
using namespace std;


class Log {

public:
  enum Level {
    TRACE = 0, DEBUG = 10, INFO = 20, WARN = 30, ERROR = 40 
  };
  
public:
  Log(const string& name) 
    : _name(name), _level(INFO), _nostream(new ostream(0)), _writeTime(true) { }

  Log(const string& name, const Level& level) 
    : _name(name), _level(level), _nostream(new ostream(0)), _writeTime(true) { }

public:
  static Log& getLog(const string& name) {
    return *(new Log(name));
  }

  static Log& getLog(const string& name, const Level& level) {
    return *(new Log(name, level));
  }
    
public:
  Level getLevel() const {
    return _level;
  }

  Log& setLevel(const Level& level) {
    _level = level;
    return *this;
  }

  static string getLevelName(const Level& level) {
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

  Log& setName(const string& name) {
    _name = name;
    return *this;
  }

  bool writeTime() const {
    return _writeTime;
  }

  Log& writeTime(const bool writeTime) {
    _writeTime = writeTime;
    return *this;
  }

  bool isActive(const Level& level) const {
    return (level >= _level);
  }

private:
  string _name;

  Level _level;

  bool _writeTime;

  ostream* const _nostream;

public:
  friend ostream& operator<<(Log& log, const Level& level);

};
  
  
ostream& operator<<(Log& log, const Log::Level& level) {
  if (log.isActive(level)) {
    time_t rawtime;
    time(&rawtime);
    char* timestr = ctime(&rawtime);
    timestr[24] = ' ';

    cout << Log::getLevelName(level) << " ";
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
  Log log = Log::getLog("foo");
  log << Log::INFO << "hello" << endl;
  log << Log::DEBUG << "hi" << endl;
  log.setLevel(Log::DEBUG);
  log << Log::WARN << "hola" << endl;
  log << Log::DEBUG << "hey" << endl;
  return EXIT_SUCCESS;
}
