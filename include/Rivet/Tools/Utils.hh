// -*- C++ -*-
#ifndef RIVET_Utils_HH
#define RIVET_Utils_HH

#include <Rivet/Math/Math.hh>
#include <Rivet/RivetBoost.hh>
#include <cctype>
#include <algorithm>
#include <cerrno>


namespace Rivet {


  /// @brief Convert any object to a string
  ///
  /// Just a convenience wrapper for the more general Boost lexical_cast
  template <typename T>
  inline std::string to_str(const T& x) {
    return lexical_cast<string>(x);
  }

  /// @brief Convert any object to a string
  ///
  /// An alias for to_str() with a more "Rivety" mixedCase name.
  template <typename T>
  inline std::string toString(const T& x) {
    return to_str(x);
  }


  /// Case-insensitive string comparison function
  ///
  /// @todo Replace with something from the Boost string library
  inline int nocase_cmp(const std::string& s1, const std::string& s2) {
    string::const_iterator it1 = s1.begin();
    string::const_iterator it2 = s2.begin();
    while ( (it1 != s1.end()) && (it2 != s2.end()) ) {
      if(::toupper(*it1) != ::toupper(*it2)) { // < Letters differ?
        // Return -1 to indicate smaller than, 1 otherwise
        return (::toupper(*it1) < ::toupper(*it2)) ? -1 : 1;
      }
      // Proceed to the next character in each string
      ++it1;
      ++it2;
    }
    size_t size1 = s1.size(), size2 = s2.size(); // Cache lengths
    // Return -1,0 or 1 according to strings' lengths
    if (size1 == size2) return 0;
    return (size1 < size2) ? -1 : 1;
  }


  /// Case-insensitive string equality function
  inline bool nocase_equals(const std::string& s1, const std::string& s2) {
    return nocase_cmp(s1, s2) == 0;
  }


  /// Convert a string to lower-case
  ///
  /// @todo Replace with something from the Boost string library
  inline std::string toLower(const std::string& s) {
    string out = s;
    transform(out.begin(), out.end(), out.begin(), (int(*)(int)) tolower);
    return out;
  }


  /// Convert a string to upper-case
  ///
  /// @todo Replace with something from the Boost string library
  inline std::string toUpper(const std::string& s) {
    string out = s;
    std::transform(out.begin(), out.end(), out.begin(), (int(*)(int)) toupper);
    return out;
  }


  /// Check whether a string @a start is found at the start of @a s
  inline bool startsWith(const std::string& s, const std::string& start) {
    if (s.length() < start.length()) return false;
    return s.substr(0, start.length()) == start;
  }


  /// Check whether a string @a end is found at the end of @a s
  inline bool endsWith(const std::string& s, const std::string& end) {
    if (s.length() < end.length()) return false;
    return s.substr(s.length() - end.length()) == end;
  }


  /// @brief Split a path string with colon delimiters.
  ///
  /// Ignores zero-length substrings. Designed for getting elements of filesystem paths, naturally.
  inline vector<std::string> pathsplit(const std::string& path) {
    const string delim = ":";
    vector<string> dirs;
    string tmppath = path;
    while (true) {
      const size_t delim_pos = tmppath.find(delim);
      if (delim_pos == string::npos) break;
      const string dir = tmppath.substr(0, delim_pos);
      if (dir.length()) dirs.push_back(dir); // Don't insert "empties"
      tmppath.replace(0, delim_pos+1, "");
    }
    if (tmppath.length()) dirs.push_back(tmppath); // Don't forget the trailing component!
    return dirs;
  }


  /// @brief Join several filesystem paths together with a delimiter character.
  ///
  /// Note that this does NOT join path elements together with a platform-portable
  /// directory delimiter, cf. the Python @c {os.path.join}!
  inline std::string pathjoin(const vector<std::string>& paths) {
    const string delim = ":";
    string rtn;
    for (vector<string>::const_iterator is = paths.begin(); is != paths.end(); ++is) {
      if (rtn.size() > 0) rtn += delim;
      rtn += *is;
    }
    return rtn;
  }


}
#endif



#ifndef CEDARSTD
#define CEDARSTD
namespace std {


  /// @name Standard library enhancements
  //@{


  /// @todo Add a "contains" function for vectors, maps and strings (returns bool rather than iterator)

  /// @todo Merge in nice ideas from mili: https://code.google.com/p/mili/


  template <typename T>
  inline void operator+=(set<T>& s1, const set<T>& s2) {
    for (typename set<T>::const_iterator s = s2.begin(); s != s2.end(); ++s) {
      s1.insert(*s);
    }
  }

  template <typename T>
  inline set<T> operator+(const set<T>& s1, const set<T>& s2) {
    set<T> rtn(s1);
    rtn += s2;
    return rtn;
  }

  template <typename T>
  inline std::string join(const set<T>& s, const std::string& sep = " ") {
    stringstream out;
    bool first = false;
    for (typename set<T>::const_iterator it = s.begin(); it != s.end(); ++it) {
      if (first) {
        first = false;
      } else {
        out << sep;
      }
      out << *it;
    }
    return out.str();
  }


  template <typename T>
  inline void operator+=(vector<T>& v1, const vector<T>& v2) {
    for (typename vector<T>::const_iterator s = v2.begin(); s != v2.end(); ++s) {
      v1.push_back(*s);
    }
  }

  template <typename T>
  inline vector<T> operator+(const vector<T>& v1, const vector<T>& v2) {
    vector<T> rtn(v1);
    rtn += v2;
    return rtn;
  }

  template <typename T>
  inline string join(const vector<T>& v, const std::string& sep = " ") {
    stringstream out;
    for (size_t i = 0; i < v.size(); ++i) {
      if (i != 0) out << sep;
      out << v[i];
    }
    return out.str();
  }

  //@}


}
#endif
