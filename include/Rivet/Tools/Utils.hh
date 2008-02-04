// -*- C++ -*-
#ifndef RIVET_Utils_HH
#define RIVET_Utils_HH

#include <Rivet/Rivet.hh>
#include <Rivet/Math/Math.hh>
#include <cctype>
#include <algorithm>
#include <cerrno>


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


  /// Split a string with single-character delimiters, ignoring zero-length 
  /// substrings. Designed for getting elements of filesystem paths, naturally.
  inline vector<string> split(string path, const string delim = ":") {
    vector<string> dirs;
    if (delim.length() != 1) {
      throw runtime_error("Rivet::split(string): delimiter must be a single character.");
    }
    while (true) {
      const size_t delim_pos = path.find(delim);
      if (delim_pos == string::npos) break;
      const string dir = path.substr(0, delim_pos);
      if (dir.length()) dirs.push_back(dir); // Don't insert "empties"
      path.replace(0, delim_pos+1, "");
    }
    if (path.length()) dirs.push_back(path); // Don't forget the trailing component!
    return dirs;
  }

  /// Get library install path
  const string getLibPath();

  /// Get data install path
  const string getDataPath();

  /// Get Rivet data install path
  const string getRivetDataPath();

  /// Get RivetGun data install path
  const string getRivetgunDataPath();


  // Return distance of closest approach from track to given (primary) vertex position.
  double get2dClosestApproach(const HepMC::GenParticle& track, const Vector3& vtx3pos);

  // Return distance of closest approach from track to given (primary) vertex position.
  double get3dClosestApproach(const HepMC::GenParticle& track, const Vector3& vtx3pos);

  /// Return 2-dimensional decay length between two vertices in transverse plane.
  double get2dDecayLength(const Vector3& vtx1, const Vector3& vtx2, const FourMomentum& jetaxis);

  /// Return 3-dimensional decay length between vertices.
  double get3dDecayLength(const Vector3& vtx1, const Vector3& vtx2, const FourMomentum& jetaxis);

}
#endif


#ifndef CEDARSTD
#define CEDARSTD
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

  template <typename T>
  inline string join(const vector<T>& v, const string& sep = " ") {
    stringstream out; 
    for (size_t i = 0; i < v.size(); ++i) {
      if (i != 0) out << sep;
      out << v[i];
    }
    return out.str();
  }

}
#endif
