#include "Rivet/AnalysisInfo.hh"
#include "yaml.h"

namespace Rivet {


  /// Ideas: 
  ///  * store as pointer on Analysis: populate only when requested
  ///  * search RIVET_DATA_PATH for <name>.info.yaml
  ///  * AnalysisInfo::make(name) returns (possibly null) pointer
  ///    - require smart pointer to delete automatically when Analysis 
  ///      goes out of scope


  /// Static factory method
  const AnalysisInfo* AnalysisInfo::make(const std::string& name) {
    /// @todo Search metadata path and read first matching file
    //if (notfound) return 0;

    /// @todo Get document

    AnalysisInfo* ai = new AnalysisInfo();
    ai->_name = "NAME";
    ai->_spiresId = "12345678";
    ai->_authors += "foo <bar@baz.tld";
    ai->_summary = "blah";
    // ai->_description;
    // ai->_runInfo;
    // ai->_experiment;
    // ai->_collider;
    // ai->_beams;
    ai->_year = "2009";
    //ai->_references;
    ai->_status = "VALIDATED";
    ai->_needsCrossSection = false;

    /// @todo Push values into read-only private members
    return ai;
  }


  string toString(const AnalysisInfo& ai) {
    /// @todo Fill in the gap...
    return "TODO";
  }

}
