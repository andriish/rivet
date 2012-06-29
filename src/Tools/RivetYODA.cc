#include "Rivet/RivetYODA.hh"
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/RivetPaths.hh"

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
    /// @todo Remove debug cerr
    cerr << "HERE2 " << aovec.size() << '\n';

    // Return value, to be populated
    RefDataMap rtn;
    foreach ( YODA::AnalysisObject * ao, aovec ) {
      Scatter2DPtr refdata( dynamic_cast<Scatter2D *>(ao) );
      if (!refdata) continue;
      string plotpath = refdata->path();
      /// @todo Remove debug cerr
      cerr << plotpath << '\n';
      rtn[plotpath] = refdata;
    }
    return rtn;
  }


}
