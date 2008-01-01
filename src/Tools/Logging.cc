
#include "Rivet/Tools/Logging.hh"
#include <ctime>

using namespace std;


namespace Rivet {

  Log::LogMap Log::existingLogs;
  Log::LevelMap Log::defaultLevels;
  Log::ColorCodes Log::colorCodes;
  string Log::endColorCode;
  bool Log::showTimestamp = false;
  bool Log::showLogLevel = true;
  bool Log::showLoggerName = true;
  bool Log::useShellColors = true;


  Log::Log(const string& name) 
    : _name(name), _level(INFO), _nostream(new ostream(0)) { }


  Log::Log(const string& name, const Level& level)
    : _name(name), _level(level), _nostream(new ostream(0)) { }


  Log& Log::getLog(const string& name) {
    if (existingLogs.find(name) == existingLogs.end()) {
      Level level = INFO;
      // Try running through all parent classes to find an existing level
      string tmpname = name;
      bool triedAllParents = false;
      while (! triedAllParents) {
        // Is there a default level?
        if (defaultLevels.find(tmpname) != defaultLevels.end()) {
          level = defaultLevels.find(tmpname)->second;
          break;
        }
        // Is there already such a logger?
        if (existingLogs.find(tmpname) != existingLogs.end()) {
          level = existingLogs.find(tmpname)->second->getLevel();
          break;
        }
        // Crop the string back to the next parent level
        size_t lastDot = tmpname.find_last_of(".");
        if (lastDot == string::npos) {
          triedAllParents = true;
        } else {
          tmpname = tmpname.substr(0, lastDot);
        }
      }
      existingLogs[name] = new Log(name, level);
    }
    return *existingLogs[name];
  }


  Log& Log::getLog(const string& name, const Level& level) {
    //cout << "Log::getLog(name, level): are we really using this? Does it make sense?" << endl;
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


  string Log::getColorCode(const Level& level) {
    if (!Log::useShellColors) return "";
    if (Log::colorCodes.empty()) {
        char* env = 0;
        /// @todo Is this "shell" (i.e. lower case) for (t)csh shells?
        /// @todo Should we actually be testing for $TERM?
        env = getenv("SHELL");
        bool invalidShell = true;

        // If there is a valid shell, try to use the appropriate codes.
        if (env) {
          //const string shell = env;
          /// @todo In a perfect world, these codes could be customised via environment / cfg file...
          /// @todo The bash codes seem to work with tcsh... for me at least (AB)
          //if (shell.find("/bash") != string::npos) {
            invalidShell = false;
            Log::colorCodes[TRACE] = "\033[0;37m";
            Log::colorCodes[DEBUG] = "\033[0;36m";
            Log::colorCodes[INFO]  = "\033[0;32m";
            Log::colorCodes[WARN]  = "\033[0;33m";
            Log::colorCodes[ERROR] = "\033[0;31m";
            Log::endColorCode      = "\033[0m";
          /// @todo Test this on csh, zsh, tcsh etc. Ncurses is too much pain to be worth trying.
          // } else if (shell.find("/tcsh") != string::npos) {
          //   invalidShell = false;
          //   Log::colorCodes[TRACE] = "";
          //   Log::colorCodes[DEBUG] = "";
          //   Log::colorCodes[INFO]  = "";
          //   Log::colorCodes[WARN]  = "";
          //   Log::colorCodes[ERROR] = "";
          //   Log::endColorCode      = "";
          // }
        }

        // If there is no such environment variable or it's unknown,
        // fill the map with empty codes.
        if (invalidShell) {
          Log::colorCodes[TRACE] = "";
          Log::colorCodes[DEBUG] = "";
          Log::colorCodes[INFO] = "";
          Log::colorCodes[WARN] = "";
          Log::colorCodes[ERROR] = "";
        }
      }

      // Return the appropriate code from the colour map.
      return colorCodes[level];
    }


  Log::Level Log::getLevelFromName(const string& level) {
    if (level == "TRACE") return TRACE;
    if (level == "DEBUG") return DEBUG;
    if (level == "INFO") return INFO;
    if (level == "WARN") return WARN;
    if (level == "ERROR") return ERROR;
    throw runtime_error("Couldn't create a log level from string '" + level + "'");
  }


  string Log::formatMessage(const Level& level, const string& message) {
    string out;
    if (Log::useShellColors) {
      out += getColorCode(level);
    }

    if (Log::showLoggerName) {
      out += getName();
      out += ": ";
    }

    if (Log::showLogLevel) {
      out += Log::getLevelName(level);
      out += " ";
    }

    if (Log::showTimestamp) {
      time_t rawtime;
      time(&rawtime);
      char* timestr = ctime(&rawtime);
      timestr[24] = ' ';
      out += timestr;
      out += " ";
    }

    if (Log::useShellColors) {
      out += endColorCode;
    }

    out += " ";
    out += message;
    
    return out;
  }


  void Log::log(const Level& level, const string& message) {
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
