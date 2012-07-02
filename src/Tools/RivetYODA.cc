#include "Rivet/RivetYODA.hh"
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/RivetPaths.hh"
#include "boost/algorithm/string/split.hpp"

using namespace std;

namespace Rivet {


  string getDatafilePath(const string& papername) {
    const string path =  findAnalysisRefFile(papername + ".aida");
    if (!path.empty()) return path;
    throw Rivet::Error("Couldn't find ref data file '" + papername + ".aida" +
                       " in $RIVET_REF_PATH, " + getRivetDataPath() + ", or .");
    return "";
  }


  RefDataMap getRefData(const string& papername) {
    // Get filename
    const string xmlfile = getDatafilePath(papername);

    YODA::Reader & reader =  ReaderAIDA::create();
    vector<YODA::AnalysisObject *> aovec;
    reader.read(xmlfile, aovec);

    // Return value, to be populated
    RefDataMap rtn;
    foreach ( YODA::AnalysisObject * ao, aovec ) {
      Scatter2DPtr refdata( dynamic_cast<Scatter2D *>(ao) );
      if (!refdata) continue;
      string plotpath = refdata->path();

      // split path at "/" and only return the last field, i.e. the histogram ID
      std::vector<string> pathvec;
      split( pathvec, plotpath, is_any_of("/"), token_compress_on );
      plotpath = pathvec.back();

      rtn[plotpath] = refdata;
    }
    return rtn;
  }


}
