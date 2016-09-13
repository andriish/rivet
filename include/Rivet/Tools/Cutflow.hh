#ifndef RIVET_Cutflow_HH
#define RIVET_Cutflow_HH

#include "Rivet/Tools/Utils.hh"

namespace Rivet {


  /// A tracker of numbers & fractions of events passing sequential cuts
  struct Cutflow {

    /// @todo Also provide single-cut filling/reading by index/name
    /// @todo Provide a Cutflows wrapper with index/name access to contained CFs and index/name all-CF filling
    /// @todo Print total and incremental acceptance percentages

    /// @brief Default constructor
    ///
    /// Does nothing! Just to allow storage in STL containers and use as a member variable without using the init list
    Cutflow() {}

    /// Proper constructor
    Cutflow(const string& cfname, const vector<string>& cutnames)
      : name(cfname), ncuts(cutnames.size()), cuts(cutnames), counts(ncuts+1, 0)
    {  }

    /// @brief Fill the pre-cut counter
    void fillinit() {
      counts[0] += 1;
    }

    /// @brief Fill the @a {icut}'th post-cut counter
    bool fill(size_t icut, bool cutresult) {
      if (cutresult) counts[icut+1] += 1;
      return cutresult;
    }

    /// @brief Fill all cut-state counters from an Ncut-element results vector
    ///
    /// This function is to be used to fill all of an event's pre- and post-cut
    /// state counters at once, including the incoming event counter. It must not be
    /// mixed with calls to the @c fill(size_t, bool) and @c fillinit() methods,
    /// or double-counting will occur.
    bool fill(const vector<bool>& cutresults) {
      if (cutresults.size() != ncuts)
        throw RangeError("Number of filled cut results needs to match the Cutflow construction");
      counts[0] += 1;
      for (size_t i = 0; i < ncuts; ++i) {
        if (cutresults[i]) counts[i+1] += 1; else break;
      }
      return all(cutresults, [](const bool x){return x;});
    }

    /// @brief Fill the N trailing post-cut counters, when supplied with an N-element results vector
    ///
    /// The @a cutresults vector represents the boolean results of the last N cuts. This function
    /// allows mixing of cut-flow filling with higher-level analyze() function escapes such as
    /// the vetoEvent directive. The initial state (state 0) is not incremented.
    bool filltail(const vector<bool>& cutresults) {
      if (cutresults.size() > ncuts)
        throw RangeError("Number of filled cut results needs to match the Cutflow construction");
      const size_t offset = counts.size() - cutresults.size();
      for (size_t i = 0; i < cutresults.size(); ++i) {
        if (cutresults[i]) counts[offset+i] += 1; else break;
      }
      return all(cutresults);
    }

    /// Create a string representation
    string str() const {
      stringstream ss;
      ss << name << " cut-flow:";
      size_t maxlen = 0;
      for (const string& t : cuts) maxlen = max(t.length(), maxlen);
      for (size_t i = 0; i <= ncuts; ++i) {
        const int pcttot = (counts[0] == 0) ? -1 : int(100*counts[i]/double(counts[0]));
        const int pctinc = (i == 0 || counts[i-1] == 0) ? -1 : int(100*counts[i]/double(counts[i-1]));
        ss << "\n" << setw(maxlen+5) << left
           << (i == 0 ? "" : "Pass "+cuts[i-1]) << "   " << right
           << setw(toString(counts[0]).length()) << toString(counts[i]) << "    "
           << setw(4) << (pcttot < 0 ? "- " : toString(pcttot)+"%") << "    "
           << setw(4) << (pctinc < 0 ? "- " : toString(pctinc)+"%");
      }
      return ss.str();
    }

    /// Print string representation to a stream
    void print(ostream& os) const {
      os << str() << flush;
    }

    string name;
    size_t ncuts;
    vector<string> cuts;
    vector<int> counts;

  };

  /// Print a Cutflow to a stream
  ostream& operator << (ostream& os, const Cutflow& cf) {
    return os << cf.str();
  }



  /// A container for several Cutflow objects, with some convenient batch access
  struct Cutflows {

    /// Do-nothing default constructor
    Cutflows() {  }

    /// Populating constructor
    Cutflows(const vector<Cutflow>& cutflows) : cfs(cutflows) {  }

    /// Append a provided Cutflow to the list
    void addCutflow(const Cutflow& cf) {
      cfs.push_back(cf);
    }

    /// Append a newly constructed Cutflow to the list
    void addCutflow(const string& cfname, const vector<string>& cutnames) {
      cfs.push_back(Cutflow(cfname, cutnames));
    }

    /// Access the @a i'th Cutflow
    Cutflow& operator [] (size_t i) { return cfs[i]; }
    /// Access the @a i'th Cutflow (const)
    const Cutflow& operator [] (size_t i) const { return cfs[i]; }

    /// Access the Cutflow whose name is @a name
    Cutflow& operator [] (const string& name) {
      for (Cutflow& cf : cfs)
        if (cf.name == name) return cf;
      throw UserError("Requested cut-flow name '" + name + "' does not exist");
    }
    /// Access the @a i'th Cutflow (const)
    const Cutflow& operator [] (const string& name) const {
      for (const Cutflow& cf : cfs)
        if (cf.name == name) return cf;
      throw UserError("Requested cut-flow name '" + name + "' does not exist");
    }

    /// Fill the pre-cuts state counter for all contained Cutflows
    void fillinit() {
      for (Cutflow& cf : cfs) cf.fillinit();
    }

    /// Create a string representation
    string str() const {
      stringstream ss;
      for (const Cutflow& cf : cfs)
        ss << cf << "\n\n";
      return ss.str();
    }

    /// Print string representation to a stream
    void print(ostream& os) const {
      os << str() << flush;
    }

    vector<Cutflow> cfs;

  };

  /// Print a Cutflows to a stream
  ostream& operator << (ostream& os, const Cutflows& cfs) {
    return os << cfs.str();
  }


}

#endif
