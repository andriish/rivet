// -*- C++ -*-
#ifndef RIVET_Utils_HH
#define RIVET_Utils_HH 1

#include <cctype>
#include <algorithm>

namespace Rivet {

  inline int nocase_cmp(const string& s1, const string& s2) {
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


  inline string toLower(const string& s) {
    string out = s;
    transform(out.begin(), out.end(), out.begin(), (int(*)(int)) tolower); 
    return out;
  }


  inline string toUpper(const string& s) {
    string out = s;
    std::transform(out.begin(), out.end(), out.begin(), (int(*)(int)) toupper); 
    return out;
  }


  template <typename Real>
  inline bool fuzzyEquals(Real a, Real b, Real tolerance = 0.001) {
    const double absavg = fabs(a + b)/2.0;
    const double absdiff = fabs(a - b);
    return (absavg == 0.0 && absdiff == 0.0) || absdiff/absavg < tolerance;
  }


  /// Split a string with single-character delimiters, ignoring zero-length 
  /// substrings. Designed for getting elements of filesystem paths, naturally.
  inline vector<string> split(const string& s, const string delim = ":") {
    string path = s;
    vector<string> dirs;
    if (delim.length() != 1) {
      throw runtime_error("Rivet::split(string): delimiter must be a single character.");
    }
    while (size_t delim_pos = path.find(delim) != string::npos) {
      string dir = path.substr(0, delim_pos);
      if (dir.length()) dirs.push_back(dir); // Don't insert "empties"
      path.replace(0, delim_pos+1, "");
    }
    if (path.length()) dirs.push_back(path); // Don't forget the trailing component!
    return dirs;
  }


  /// Get Rivet data install path
  const string getInstalledDataPath();

  /// Get Rivet library install path
  const string getInstalledLibPath();

}


namespace std {

  template <typename T>
  inline void operator+=(set<T>& s1, const set<T>& s2) {
    for (typename set<T>::const_iterator s = s2.begin(); s != s2.end(); ++s) {
      s1.insert(*s);
    }
  }

  template <typename T>
  inline void operator+=(vector<T>& s1, const vector<T>& s2) {
    for (typename vector<T>::const_iterator s = s2.begin(); s != s2.end(); ++s) {
      s1.push_back(*s);
    }
  }

}

#endif
