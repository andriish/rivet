// $Id: $
#ifndef RIVET_LOGGING_H 
#define RIVET_LOGGING_H 1

#include <string>
#include <iostream>

namespace Rivet {
  
  class Log {
  public:
    enum Level {
      TRACE = 0, DEBUG = 10, INFO = 20, WARN = 30, ERROR = 40 
    };

  private:
    Log(const std::string& name);

    Log(const std::string& name, const Level& level);

  public:
    static Log& getLog(const std::string& name);

    static Log& getLog(const std::string& name, const Level& level);

  public:
    inline Level getLevel() const {
      return _level;
    }

    inline Log& setLevel(const Level& level) {
      _level = level;
      return *this;
    }

    static std::string getLevelName(const Level& level);

    inline std::string getName() const {
      return _name;
    }

    inline Log& setName(const std::string& name) {
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

    inline void trace(const std::string& message) { log(TRACE, message); }

    inline void debug(const std::string& message) { log(DEBUG, message); }
    
    inline void info(const std::string& message) { log(INFO, message); }

    inline void warn(const std::string& message) { log(WARN, message); }

    inline void error(const std::string& message) { log(ERROR, message); }

  private:
    std::string _name;
    
    Level _level;
    
    bool _writeTime;
    
  protected:
    //public:

    void log(const Level& level, const std::string& message);

    std::string formatMessage(const Level& level, const std::string& message);

  public:

    /// @todo Hide this...
    std::ostream* const _nostream;

    friend std::ostream& operator<<(Log& log, const Log::Level& level);

  };
  

  std::ostream& operator<<(Log& log, const Log::Level& level);
  
}


#endif // RIVET_LOGGING_H
