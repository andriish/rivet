// $Id: $
#ifndef RIVET_LOGGING_H 
#define RIVET_LOGGING_H 1

#include "Rivet/Rivet.hh"


namespace Rivet {

  class Log {
  public:

    /// Log priority levels.
    enum Level {
      TRACE = 0, DEBUG = 10, INFO = 20, WARN = 30, ERROR = 40 
    };

    /// Typedef for a collection of named logs.
    typedef map<const string, Log*> LogMap;

    /// Typedef for a collection of named log levels.
    typedef map<const string, Level> LevelMap;

  private:
    /// A static map of existing logs: we don't make more loggers than necessary.
    static LogMap existingLogs;

    /// A static map of default log levels.
    static LevelMap defaultLevels;

    /// Show timestamp?
    static bool showTimestamp;

    /// Show log level?
    static bool showLogLevel;

    /// Show logger name?
    static bool showLoggerName;

  public:
    /// Set the default log levels
    inline static void setDefaultLevels(const LevelMap& logLevels) {
      defaultLevels = logLevels;
    }

    inline static void setShowTimestamp(const bool showTime=true) {
      showTimestamp = showTime;
    }

    inline static void setShowLevel(const bool showLevel=true) {
      showLogLevel = showLevel;
    }

    inline static void setShowLoggerName(const bool showName=true) {
      showLoggerName = showName;
    }


  protected:
    /// @name Hidden constructors etc.
    //@{
    /// Constructor 1
    Log(const string& name);

    /// Constructor 2
    Log(const string& name, const Level& level);

    /// Copy constructor
    //Log(const Log&);

    /// Copy assignment operator
    //Log& operator=(const Log&);
    //@}

  public:
    /// Get a logger with the given name. The level will be taken from the 
    /// "requestedLevels" static map or will be INFO by default.
    static Log& getLog(const string& name);

    /// Get a logger with the given name and level.
    static Log& getLog(const string& name, const Level& level);

  public:
    /// Get the priority level of this logger.
    inline Level getLevel() const {
      return _level;
    }

    /// Set the priority level of this logger.
    inline Log& setLevel(const Level& level) {
      _level = level;
      return *this;
    }

    /// Get a log level enum from a string.
    static Level getLevelFromName(const string& level);

    /// Get the string representation of a log level.
    static string getLevelName(const Level& level);

    /// Get the name of this logger.
    inline string getName() const {
      return _name;
    }

    /// Set the name of this logger.
    inline Log& setName(const string& name) {
      _name = name;
      return *this;
    }

    /// Will this log level produce output on this logger at the moment?
    inline bool isActive(const Level& level) const {
      return (level >= _level);
    }

    /// @name Explicit log methods
    //@{
    inline void trace(const string& message) { log(TRACE, message); }

    inline void debug(const string& message) { log(DEBUG, message); }
    
    inline void info(const string& message) { log(INFO, message); }

    inline void warn(const string& message) { log(WARN, message); }

    inline void error(const string& message) { log(ERROR, message); }
    //@}

  private:
    /// This logger's name
    string _name;
    
    /// Threshold level for this logger.
    Level _level;
    
  protected:
    /// Write a message at a particular level.
    void log(const Level& level, const string& message);

    /// Turn a message string into the current log format.
    string formatMessage(const Level& level, const string& message);

  public:

    /// A null output stream, used for piping discarded output to nowhere.
    /// @todo Hide this...
    ostream* const _nostream;

    /// The streaming operator can use Log's internals.
    friend ostream& operator<<(Log& log, const Log::Level& level);

  };
  
  /// Streaming output to a logger must have a Log::Level as its first argument.
  ostream& operator<<(Log& log, const Log::Level& level);
  
}


#endif // RIVET_LOGGING_H
