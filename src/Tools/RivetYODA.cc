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


  map<string, YODA::AnalysisObjectPtr> getRefData(const string& papername) {
    const string datafile = getDatafilePath(papername);

    // Make an appropriate data file reader and read the data objects
    /// @todo Remove AIDA support some day...
    YODA::Reader& reader = (datafile.find(".yoda") != string::npos) ?   \
      YODA::ReaderYODA::create() : YODA::ReaderAIDA::create();
    vector<YODA::AnalysisObject *> aovec;
    reader.read(datafile, aovec);

    // Return value, to be populated
    map<string, YODA::AnalysisObjectPtr> rtn;
    foreach ( YODA::AnalysisObject* ao, aovec ) {
      YODA::AnalysisObjectPtr refdata(ao);
      if (!refdata) continue;
      const string plotpath = refdata->path();
      // Split path at "/" and only return the last field, i.e. the histogram ID
      const size_t slashpos = plotpath.rfind("/");
      const string plotname = (slashpos+1 < plotpath.size()) ? plotpath.substr(slashpos+1) : "";
      rtn[plotname] = refdata;
    }
    return rtn;
  }
}

namespace {

  using Rivet::Fill;
  using Rivet::Fills;

  template <class T>
  void fillAllPersistent(const vector<typename T::Ptr> & persistent, 
                         typename T::FillType x, double w, 
                         const vector<double> & weight) {
    for ( size_t m = 0; m < persistent.size(); ++m ) {
      persistent[m]->fill( x, w * weight[m] );
    }
  }

    template <class T>
    double get_window_size(const typename T::Ptr & histo,
                           typename T::BinType x) {
        // the bin index we fall in
        const auto binidx = histo->binIndexAt(x);
        // gaps, overflow, underflow don't contribute
        if ( binidx == -1 ) 
            return 0;


        const auto & b = histo->bin(binidx);

        // if we don't have a valid neighbouring bin,
        // we use infinite width
        typename T::Bin b1(-1.0/0.0, 1.0/0.0);

        // points in the top half compare to the upper neighbour
        if ( x > b.xMid() ) {
            int nextidx = binidx + 1;
            if ( nextidx < histo->bins().size() )
                b1 = histo->bin(nextidx);
        }
        else { // compare to the lower neighbour
            int nextidx = binidx - 1;
            if ( nextidx >= 0 )
                b1 = histo->bin(nextidx);
        }
        // the factor 2 is arbitrary, could poss. be smaller
        return min( b.width(), b1.width() ) / 2.0;
    }

    template <class T>
     typename T::BinType 
     fillT2binT(typename T::FillType a) {
       return a;
    }

    template <>
    YODA::Profile1D::BinType 
    fillT2binT<YODA::Profile1D>(YODA::Profile1D::FillType a) {
      return get<0>(a);
    }

    template <>
    YODA::Profile2D::BinType 
    fillT2binT<YODA::Profile2D>(YODA::Profile2D::FillType a) {
      return { get<0>(a), get<1>(a) };
    }



  template <class T>
  void commit(typename T::Ptr & persistent,
              const vector< vector<Fill<T> > > & tuple,
              const vector<double> weights /* generator weight over subevents */) {
    for ( const auto & x : tuple ) {
      double maxwindow = 0.0;
      for ( const auto & xi : x ) {
        // check for NOFILL here
        double window = get_window_size<T>(persistent, fillT2binT<T>(xi.first));
        if ( window > maxwindow )
          maxwindow = window;
      }
      const double wsize = maxwindow;
      // all windows have same size

      set<double> edgeset;
      // bin edges need to be in here!
      for ( const auto & xi : x ) {
        edgeset.insert(fillT2binT<T>(xi.first) - wsize);
        edgeset.insert(fillT2binT<T>(xi.first) + wsize);
      }

      vector< std::tuple<double,double,double> > hfill;
      double sumf = 0.0;
      auto edgit = edgeset.begin();
      double ehi = *edgit;
      while ( ++edgit != edgeset.end() ) {
        double elo = ehi;
        ehi = *edgit;
        double sumw = 0.0;
        bool gap = true; // Check for gaps between the sub-windows.
        for ( int i = 0; i < x.size(); ++i  ) {
          // check equals comparisons here!
          if ( fillT2binT<T>(x[i].first) + wsize >= ehi 
               && 
               fillT2binT<T>(x[i].first) - wsize <= elo ) {
            sumw += x[i].second * weights[i];
            gap = false;
          }
        }
        if ( gap ) continue;
        hfill.push_back(make_tuple((ehi + elo)/2.0, sumw, ehi - elo));
        sumf += ehi - elo;
      }
      for ( auto f : hfill )
        persistent->fill( get<0>(f), get<1>(f), get<2>(f)/sumf );
        // Note the scaling to one single fill
    }
  }

    template <class T>
    double distance(T a, T b) {
      return abs(a - b);
    }

    template <>
    double distance<tuple<double,double> >(tuple<double,double> a, tuple<double,double> b) {
      return Rivet::sqr(get<0>(a) - get<0>(b)) + Rivet::sqr(get<1>(a) - get<1>(b));
    }



}




/// fills is a vector of sub-event with an ordered set of x-values of
/// the fills in each sub-event. NOFILL should be an "impossible"
/// value for this histogram. Returns a vector of sub-events with
/// an ordered vector of fills (including NOFILLs) for each sub-event.
template <class T>
vector< vector<Fill<T> > > 
match_fills(const vector< Fills<T> > & fills, const Fill<T> & NOFILL) 
{
  vector< vector<Fill<T> > > matched;
  // First just copy subevents into vectors and find the longest vector.
  unsigned int maxfill = 0; // length of biggest vector
  int imax = 0; // index position of biggest vector
  for ( const auto & subev : fills ) {
    if ( subev.size() > maxfill ) {
      maxfill = subev.size();
      imax = matched.size();
    }
    matched.push_back(vector<Fill<T> >(subev.begin(), subev.end()));
  }
  // Now, go through all subevents with missing fills.
  const vector<Fill<T>> & full = matched[imax]; // the longest one
  for ( auto & subev : matched ) {
    if ( subev.size() == maxfill ) continue;

    // Add NOFILLs to the end;
    while ( subev.size() < maxfill ) subev.push_back(NOFILL);

    // Iterate from the back and shift all fill values backwards by
    // swapping with NOFILLs so that they better match the full
    // subevent.
    for ( int i = maxfill - 1; i >= 0; --i ) {
      if ( subev[i] == NOFILL ) continue;
      int j = i;
      while ( j + 1 < maxfill && subev[j + 1] == NOFILL &&
              distance(subev[j].first, full[j].first) > distance(subev[j].first, full[j + 1].first) ) 
      {
            swap(subev[j], subev[j + 1]);
            ++j;
      }
    }
  }
  // transpose
  vector<vector<Fill<T>>> result(maxfill,vector<Fill<T>>(matched.size()));
  for (size_t i = 0; i < matched.size(); ++i)
      for (size_t j = 0; j < maxfill; ++j)
          result.at(j).at(i) = matched.at(i).at(j);
  return result;
}






namespace Rivet {

  template <class T>
  void Wrapper<T>::pushToPersistent(const vector<vector<double> >& weight) {

      // have we had subevents at all?
      const bool have_subevents = _evgroup.size() > 1;
      if ( ! have_subevents ) {
        assert( _evgroup.size() == 1 && weight.size() == 1 );
        // simple replay of all tuple entries
        // each recorded fill is inserted into all persistent weightname histos
        for ( const auto & f : _evgroup[0]->fills() ) {
          fillAllPersistent<T>( _persistent, f.first, f.second, weight[0] );
        }
      } else {
        assert( _evgroup.size() == weight.size() );

        // All the fills across subevents
        // each item in allFills is a subevent
        vector<Fills<T>> allFills;
        for ( const auto & ev : _evgroup )
          allFills.push_back( ev->fills() );

        vector< vector<Fill<T> > > 
          linedUpXs = match_fills<T>(allFills, {typename T::FillType(), 0.0});

        // TODO check if weight is transposed here!
        for ( size_t m = 0; m < _persistent.size(); ++m ) {
          commit<T>( _persistent[m], linedUpXs, weight.at(m) );
        }
      }
      _evgroup.clear();
      _active.reset();
  }

  template <>
  void Wrapper<YODA::Counter>::pushToPersistent(const vector<vector<double> >& weight) {
    for ( size_t m = 0; m < _persistent.size(); ++m ) {
      for ( size_t n = 0; n < _evgroup.size(); ++n ) {
        for ( const auto & f : _evgroup[n]->fills() ) {
          _persistent[m]->fill( f.second * weight[n][m] );
        }
      }
    }
  }

  // explicitly instantiate all wrappers

  template class Wrapper<YODA::Histo1D>;
//  template class Wrapper<YODA::Histo2D>;
  template class Wrapper<YODA::Profile1D>;
//  template class Wrapper<YODA::Profile2D>;
  template class Wrapper<YODA::Counter>;

}
