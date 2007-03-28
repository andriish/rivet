// $Id: $
#ifndef RIVET_LOGGING_H 
#define RIVET_LOGGING_H 1

#include "Rivet/Rivet.hh"


namespace Rivet {
  
  class Log {
  public:
    enum Level {
      TRACE = 0, DEBUG = 10, INFO = 20, WARN = 30, ERROR = 40 
    };

  private:
    Log(const string& name);

    Log(const string& name, const Level& level);

  public:
    static Log& getLog(const string& name);

    static Log& getLog(const string& name, const Level& level);

  public:
    inline Level getLevel() const {
      return _level;
    }

    inline Log& setLevel(const Level& level) {
      _level = level;
      return *this;
    }

    static Level getLevelFromName(const string& level);

    static string getLevelName(const Level& level);

    inline string getName() const {
      return _name;
    }

    inline Log& setName(const string& name) {
      _name = name;
      return *this;
    }

    inline bool writeTime() const {
      return _writeTime;
    }

    inline Log& writeTime(const bool writeTime) {
      _writeTime = writeTime;
      return *this;
    }

    inline bool isActive(const Level& level) const {
      return (level >= _level);
    }

    inline void trace(const string& message) { log(TRACE, message); }

    inline void debug(const string& message) { log(DEBUG, message); }
    
    inline void info(const string& message) { log(INFO, message); }

    inline void warn(const string& message) { log(WARN, message); }

    inline void error(const string& message) { log(ERROR, message); }

  private:
    string _name;
    
    Level _level;
    
    bool _writeTime;
    
  protected:
    //public:

    void log(const Level& level, const string& message);

    string formatMessage(const Level& level, const string& message);

  public:

    /// @todo Hide this...
    ostream* const _nostream;

    friend ostream& operator<<(Log& log, const Log::Level& level);

  };
  

  ostream& operator<<(Log& log, const Log::Level& level);
  
}


#endif // RIVET_LOGGING_H
