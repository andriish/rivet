// -*- C++ -*-
#ifndef RIVET_Utils_HH
#define RIVET_Utils_HH

#include "Rivet/Tools/RivetSTL.hh"
#include "Rivet/Tools/PrettyPrint.hh"
#include <sstream>
#include <cctype>
#include <algorithm>
#include <cerrno>

namespace Rivet {


  /// @name String utils
  //@{

  struct bad_lexical_cast : public std::runtime_error {
    bad_lexical_cast(const std::string& what) : std::runtime_error(what) {}
  };

  /// @brief Convert between any types via stringstream
  template<typename T, typename U>
  T lexical_cast(const U& in) {
    try {
      std::stringstream ss;
      ss << in;
      T out;
      ss >> out;
      return out;
    } catch (const std::exception& e) {
      throw bad_lexical_cast(e.what());
    }
  }

  /// @brief Convert any object to a string
  ///
  /// Just a convenience wrapper for the more general Boost lexical_cast
  template <typename T>
  inline string to_str(const T& x) {
    return lexical_cast<string>(x);
  }

  /// @brief Convert any object to a string
  ///
  /// An alias for to_str() with a more "Rivety" mixedCase name.
  template <typename T>
  inline string toString(const T& x) {
    return to_str(x);
  }

  /// Replace the first instance of patt with repl
  inline string& replace_first(string& str, const string& patt, const string& repl) {
    if (!contains(str, patt)) return str; //< contains from RivetSTL
    str.replace(str.find(patt), patt.size(), repl);
    return str;
  }

  /// @brief Replace all instances of patt with repl
  ///
  /// @note Finding is interleaved with replacement, so the second search happens after
  /// first replacement, etc. This could lead to infinite loops and other counterintuitive
  /// behaviours if not careful.
  inline string& replace_all(string& str, const string& patt, const string& repl) {
    if (!contains(str, patt)) return str; //< contains from RivetSTL
    while (true) {
      string::size_type it = str.find(patt);
      if (it == string::npos) break;
      str.replace(it, patt.size(), repl);
    }
    return str;
  }


  /// Case-insensitive string comparison function
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


  /// Case-insensitive string equality function
  inline bool nocase_equals(const string& s1, const string& s2) {
    return nocase_cmp(s1, s2) == 0;
  }


  /// Convert a string to lower-case
  inline string toLower(const string& s) {
    string out = s;
    std::transform(out.begin(), out.end(), out.begin(), (int(*)(int)) tolower);
    return out;
  }


  /// Convert a string to upper-case
  inline string toUpper(const string& s) {
    string out = s;
    std::transform(out.begin(), out.end(), out.begin(), (int(*)(int)) toupper);
    return out;
  }


  /// Check whether a string @a start is found at the start of @a s
  inline bool startsWith(const string& s, const string& start) {
    if (s.length() < start.length()) return false;
    return s.substr(0, start.length()) == start;
  }


  /// Check whether a string @a end is found at the end of @a s
  inline bool endsWith(const string& s, const string& end) {
    if (s.length() < end.length()) return false;
    return s.substr(s.length() - end.length()) == end;
  }


  /// Make a string containing the string representations of each item in v, separated by sep
  template <typename T>
  inline string join(const vector<T>& v, const string& sep=" ") {
    string rtn;
    for (size_t i = 0; i < v.size(); ++i) {
      if (i != 0) rtn += sep;
      rtn += to_str(v[i]);
    }
    return rtn;
  }

  /// Make a string containing the string representations of each item in s, separated by sep
  template <typename T>
  inline string join(const set<T>& s, const string& sep=" ") {
    string rtn;
    for (const T& x : s) {
      if (rtn.size() > 0) rtn += sep;
      rtn += to_str(x);
    }
    return rtn;
  }

  //@}


  /// @name Path utils
  //@{

  /// @brief Split a path string with colon delimiters
  ///
  /// Ignores zero-length substrings. Designed for getting elements of filesystem paths, naturally.
  inline vector<string> pathsplit(const string& path) {
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


  /// @brief Join several filesystem paths together with the standard ':' delimiter
  ///
  /// Note that this does NOT join path elements together with a platform-portable
  /// directory delimiter, cf. the Python @c {os.path.join}!
  inline string pathjoin(const vector<string>& paths) {
    return join(paths, ":");
  }

  //@}


  /// @name Container utils
  //@{

  /// Return true if f(x) is true for any x in container c, otherwise false.
  template <typename CONTAINER, typename FN>
  inline bool any(const CONTAINER& c, const FN& f) {
    for (const auto& x : c)
      if (f(x)) return true;
    return false;
  }

  /// Return true if @a f(x) is true for all @c x in container @a c, otherwise false.
  template <typename CONTAINER, typename FN>
  inline bool all(const CONTAINER& c, const FN& f) {
    for (const auto& x : c)
      if (!f(x)) return false;
    return true;
  }

  /// Generic sum function, adding @c x for all @c x in container @a c, starting with @a start
  template <typename CONTAINER, typename T>
  inline T sum(const CONTAINER& c, const T& start=T()) {
    T rtn = start;
    for (const auto& x : c) rtn += x;
    return rtn;
  }

  /// Generic sum function, adding @a fn(@c x) for all @c x in container @a c, starting with @a start
  template <typename CONTAINER, typename FN, typename T>
  inline T sum(const CONTAINER& c, const FN& f, const T& start=T()) {
    T rtn = start;
    for (const auto& x : c)
      rtn += f(x);
    return rtn;
  }

  /// A single-container-arg version of std::transform, aka @c map
  template <typename C1, typename C2, typename FN>
  inline const C2& transform(const C1& in, C2& out, const FN& f) {
    out.clear(); out.resize(in.size());
    std::transform(in.begin(), in.end(), out.begin(), f);
    return out;
  }

  /// A single-container-arg version of std::accumulate, aka @c reduce
  template <typename C1, typename T, typename FN>
  inline T accumulate(const C1& in, const T& init, const FN& f) {
    const T rtn = accumulate(in.begin(), in.end(), init, f);
    return rtn;
  }

  /// @brief Slice of the container elements cf. Python's [i:j] syntax
  ///
  /// The element at the @j index is not included in the returned container.
  /// @a i and @a j can be negative, treated as backward offsets from the end of the container.
  template <typename CONTAINER>
  inline CONTAINER slice(const CONTAINER& c, int i, int j) {
    CONTAINER rtn;
    const size_t off1 = (i >= 0) ? i : c.size() + i;
    const size_t off2 = (j >= 0) ? j : c.size() + j;
    if (off1 > c.size() || off2 > c.size()) throw RangeError("Attempting to slice beyond requested offsets");
    if (off2 < off1) throw RangeError("Requested offsets in invalid order");
    rtn.resize(off2 - off1);
    std::copy(c.begin()+off1, c.begin()+off2, rtn.begin());
    return rtn;
  }

  /// @brief Tail slice of the container elements cf. Python's [i:] syntax
  ///
  /// Single-index specialisation of @c slice(c, i, j)
  template <typename CONTAINER>
  inline CONTAINER slice(const CONTAINER& c, int i) {
    return slice(c, i, c.size());
  }

  /// @brief Head slice of the @a n first container elements
  template <typename CONTAINER>
  inline CONTAINER head(const CONTAINER& c, size_t n) {
    // if (n > c.size()) throw RangeError("Requested head longer than container");
    const size_t m = min(n, c.size());
    return slice(c, 0, m);
  }

  /// @brief Tail slice of the @a n last container elements
  template <typename CONTAINER>
  inline CONTAINER tail(const CONTAINER& c, size_t n) {
    // if (n > c.size()) throw RangeError("Requested tail longer than container");
    const size_t m = min(n, c.size());
    return slice(c, c.size()-m);
  }

  /// Find the minimum value in the vector
  double min(const vector<double>& in) {
    return *min_element(in.begin(), in.end());
    // static const auto FDBLMIN = [](double a, double b){ return std::min(a,b); };
    // const double rtn = accumulate(in.begin(), in.end(), DBL_MAX, FDBLMIN);
    // return rtn;
  }

  /// Find the maximum value in the vector
  inline double max(const vector<double>& in) {
    return *max_element(in.begin(), in.end());
    // static const auto FDBLMAX = [](double a, double b){ return std::max(a,b); };
    // const double rtn = accumulate(in.begin(), in.end(), -DBL_MAX, FDBLMAX);
    // return rtn;
  }

  /// Find the minimum and maximum values in the vector
  inline pair<double,double> minmax(const vector<double>& in) {
    const auto minmax = minmax_element(in.begin(), in.end());
    return make_pair(*minmax.first, *minmax.second);
    // static const auto FDBLMINMAX = [](double a, double b){ return std::minmax(a,b); };
    // static const auto DBL_MAXMIN = make_pair(DBL_MAX,-DBL_MAX);
    // const auto rtn = accumulate(in.begin(), in.end(), DBL_MAXMIN, FDBLMINMAX);
    // return rtn;
  }

  /// Find the minimum value in the vector
  inline int min(const vector<int>& in) {
    return *min_element(in.begin(), in.end());
    // static const auto FINTMIN = [](int a, int b){ return std::min(a,b); };
    // const int rtn = accumulate(in.begin(), in.end(), INT_MAX, FINTMIN);
    // return rtn;
  }

  /// Find the maximum value in the vector
  inline int max(const vector<int>& in) {
    return *max_element(in.begin(), in.end());
    // static const auto FINTMAX = [](int a, int b){ return std::max(a,b); };
    // const int rtn = accumulate(in.begin(), in.end(), -INT_MAX, FINTMAX);
    // return rtn;
  }

  /// Find the minimum and maximum values in the vector
  inline pair<int,int> minmax(const vector<int>& in) {
    const auto minmax = minmax_element(in.begin(), in.end());
    return make_pair(*minmax.first, *minmax.second);
    // static const auto FINTMINMAX = [](int a, int b){ return std::minmax(a,b); };
    // static const auto INT_MAXMIN = make_pair(INT_MAX,-INT_MAX);
    // const auto rtn = accumulate(in.begin(), in.end(), INT_MAXMIN, FINTMINMAX);
    // return rtn;
  }

  //@}


}

#endif
