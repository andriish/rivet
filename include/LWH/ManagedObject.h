// -*- C++ -*-
#ifndef LWH_ManagedObject_H
#define LWH_ManagedObject_H
//
// This is the declaration of the ManagedObject class.
//

#include "AIManagedObject.h"
#include <iostream>

#ifdef HAVE_ROOT
  #include "TFile.h"
#endif

namespace LWH {

using namespace AIDA;

/**
 * The creator of trees.
 */
class ManagedObject: public IManagedObject {

public:

  /// Destructor.
  virtual ~ManagedObject() {}

  /**
   * Encode sensitive characters as XML entities.
   */
  static std::string encodeForXML(const std::string& in) {
    std::string out = in;
    typedef std::pair<std::string, std::string> CharsToEntities;
    std::vector<CharsToEntities> cs2es;
    // NB. Ampersand replacement must come first!
    /// @todo Commented out ampersands until I've worked out how to avoid the infinite loop...
    // Maybe the searches should only apply to the part of the string that 
    // hasn't already been replaced...
    //cs2es.push_back(std::make_pair("&", "&amp;"));
    cs2es.push_back(std::make_pair("<", "&lt;"));
    cs2es.push_back(std::make_pair(">", "&gt;"));
    // Also quotation marks?
    for (std::vector<CharsToEntities>::const_iterator c2e = cs2es.begin(); c2e != cs2es.end(); ++c2e) {
      //std::cerr << "Found " << c2e->first << "? Pos = " << out.find(c2e->first) << std::endl;
      while (out.find(c2e->first) != std::string::npos) {
        const size_t pos = out.find(c2e->first);
        //std::cerr << "Found " << c2e->first << " at pos " << pos << " in '" << out << "'" << std::endl;
        out.replace(pos, 1, c2e->second);
      }
    }
    return out;
  }


  /**
   * Write out the object to the given stream in XML format.
   */
  virtual bool writeXML(std::ostream & os,
			std::string path, std::string name) = 0;

  /**
   * Write out the object to the given stream in simple table format.
   */
  virtual bool writeFLAT(std::ostream & os,
			 std::string path, std::string name) = 0;



#ifdef HAVE_ROOT
  /**
   * Write out the object to the given TFile in Root format.
   */
  virtual bool writeROOT(TFile* file,
			 std::string path, std::string name) = 0;
#endif


};

}

#endif /* LWH_ManagedObject_H */
