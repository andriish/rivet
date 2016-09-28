#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include "Rivet/Tools/RivetPaths.hh"
#include "YODA/ReaderYODA.h"
#include "YODA/ReaderAIDA.h"

using namespace std;

namespace Rivet {


  string getDatafilePath(const string& papername) {
    /// Try to find YODA otherwise fall back to try AIDA
    const string path1 = findAnalysisRefFile(papername + ".yoda");
    if (!path1.empty()) return path1;
    const string path2 = findAnalysisRefFile(papername + ".aida");
    if (!path2.empty()) return path2;
    throw Rivet::Error("Couldn't find ref data file '" + papername + ".yoda/aida" +
                       " in $RIVET_REF_PATH, '" + getRivetDataPath() + "', or '.'");
  }


  map<string, AnalysisObjectPtr> getRefData(const string& papername) {
    const string datafile = getDatafilePath(papername);

    // Make an appropriate data file reader and read the data objects
    /// @todo Remove AIDA support some day...
    YODA::Reader& reader = (datafile.find(".yoda") != string::npos) ?   \
      YODA::ReaderYODA::create() : YODA::ReaderAIDA::create();
    vector<YODA::AnalysisObject *> aovec;
    reader.read(datafile, aovec);

    // Return value, to be populated
    map<string, AnalysisObjectPtr> rtn;
    foreach ( YODA::AnalysisObject* ao, aovec ) {
      AnalysisObjectPtr refdata(ao);
      if (!refdata) continue;
      const string plotpath = refdata->path();
      // Split path at "/" and only return the last field, i.e. the histogram ID
      const size_t slashpos = plotpath.rfind("/");
      const string plotname = (slashpos+1 < plotpath.size()) ? plotpath.substr(slashpos+1) : "";
      rtn[plotname] = refdata;
    }
    return rtn;
  }


  void Histo1DPtr::pushToPersistent(const vector<vector<double> >& weight) {
      // loop over event weights
      for (size_t m = 0; m < _persistent.size(); ++m) {

          // this is the initial event---it always exists
          YODA::Histo1DPtr sum = boost::make_shared<YODA::Histo1D>(_evgroup[0]->clone());
          sum->scaleW( weight[0][m] );

          // loop over additional subevents (maybe there aren't any)
          // note loop starts at 1
          for (size_t n = 1; n < _evgroup.size(); ++n) {
              YODA::Histo1DPtr tmp = boost::make_shared<YODA::Histo1D>(_evgroup[n]->clone());
              tmp->scaleW( weight[n][m] );
              *sum += *tmp;
          }

          // sum typically only has one bin filled
          // do we really need to fill bin.size() times?
          bool filled_something = false;

          /// @todo
          /// this is not correct!
          /// we need to trick the histogram into thinking it was
          /// filled exactly onen time, even if subevents fall into
          /// different bins of the histogram
          /// how 
          foreach(const YODA::HistoBin1D& b, sum->bins()) {
              if ( b.effNumEntries() != 0 && b.sumW() != 0 ) {
                  // @todo number of fills will be wrong
                  // needs to be xMid. xMean is not valid if bin contains neg.weights
                  _persistent[m]->fill(b.xMid(), b.sumW());
                  filled_something = true;
              }
          }

          // @todo
          // overflow

      }

      _evgroup.clear();
      _active.reset();
  }

  void Histo2DPtr::pushToPersistent(const vector<vector<double> >& weight) {
      /// @todo

      return;
  }

  void Profile1DPtr::pushToPersistent(const vector<vector<double> >& weight) {
      /// @todo

      return;
  }

  void Profile2DPtr::pushToPersistent(const vector<vector<double> >& weight) {
      /// @todo

      return;
  }

  void Counter::pushToPersistent(const vector<vector<double> >& weight) {
      /// @todo

      return;
  }
}
