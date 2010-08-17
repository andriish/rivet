// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/AnalysisInfo.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "LWH/AIManagedObject.h"
using namespace AIDA;

namespace Rivet {


  namespace {
    string makeAxisCode(const size_t datasetId, const size_t xAxisId, const size_t yAxisId) {
      stringstream axisCode;
      axisCode << "d";
      if (datasetId < 10) axisCode << 0;
      axisCode << datasetId;
      axisCode << "-x";
      if (xAxisId < 10) axisCode << 0;
      axisCode << xAxisId;
      axisCode << "-y";
      if (yAxisId < 10) axisCode << 0;
      axisCode << yAxisId;
      return axisCode.str();
    }
  }


  ////////////////////////


  Analysis::Analysis(const string& name)
    : _crossSection(-1.0),
      _gotCrossSection(false),
      _needsCrossSection(false),
      _analysishandler(NULL),
      _madeHistoDir(false)
  {
    ProjectionApplier::_allowProjReg = false;
    _defaultname = name;
    AnalysisInfo* ai = AnalysisInfo::make(name);
    assert(ai != 0);
    _info.reset(ai);
    assert(_info.get() != 0);
    //setBeams(ANY, ANY);
  }


  Analysis::~Analysis()
  {  }


  IAnalysisFactory& Analysis::analysisFactory() {
    return handler().analysisFactory();
  }


  ITree& Analysis::tree() {
    return handler().tree();
  }


  IHistogramFactory& Analysis::histogramFactory() {
    return handler().histogramFactory();
  }


  IDataPointSetFactory& Analysis::datapointsetFactory() {
    return handler().datapointsetFactory();
  }


  double Analysis::sqrtS() const {
    return handler().sqrtS();
  }

  const ParticlePair& Analysis::beams() const {
    return handler().beams();
  }

  const PdgIdPair Analysis::beamIds() const {
    return handler().beamIds();
  }


  const string Analysis::histoDir() const {
    /// @todo This doesn't change: calc and cache at Analysis construction!
    string path = "/" + name();
    if (handler().runName().length() > 0) {
      path = "/" + handler().runName() + path;
    }
    while (find_first(path, "//")) {
      replace_all(path, "//", "/");
    }
    return path;
  }


  const string Analysis::histoPath(const string& hname) const {
    const string path = histoDir() + "/" + hname;
    return path;
  }


  Log& Analysis::getLog() const {
    string logname = "Rivet.Analysis." + name();
    return Log::getLog(logname);
  }


  size_t Analysis::numEvents() const {
    return handler().numEvents();
  }


  double Analysis::sumOfWeights() const {
    return handler().sumOfWeights();
  }


  ////////////////////////////////////////////////////////////
  // Metadata

  const AnalysisInfo& Analysis::info() const {
    assert(_info.get() != 0);
    return *_info;
  }

  string Analysis::name() const {
    if (_info && !_info->name().empty()) return _info->name();
    return _defaultname;
  }

  string Analysis::spiresId() const {
    if (!_info) return "NONE";
    return _info->spiresId();
  }

  vector<string> Analysis::authors() const {
    if (!_info) return std::vector<std::string>();
    return _info->authors();
  }

  string Analysis::summary() const {
    if (!_info) return "NONE";
    return _info->summary();
  }

  string Analysis::description() const {
    if (!_info) return "NONE";
    return _info->description();
  }

  string Analysis::runInfo() const {
    if (!_info) return "NONE";
    return _info->runInfo();
  }

  const std::vector<std::pair<double,double> >& Analysis::energies() const {
    return info().energies();
  }

  string Analysis::experiment() const {
    if (!_info) return "NONE";
    return _info->experiment();
  }

  string Analysis::collider() const {
    if (!_info) return "NONE";
    return _info->collider();
  }

  string Analysis::year() const {
    if (!_info) return "NONE";
    return _info->year();
  }

  vector<string> Analysis::references() const {
    if (!_info) return vector<string>();
    return _info->references();
  }

  string Analysis::bibKey() const {
    if (!_info) return "";
    return _info->bibKey();
  }

  string Analysis::bibTeX() const {
    if (!_info) return "";
    return _info->bibTeX();
  }

  string Analysis::status() const {
    if (!_info) return "UNVALIDATED";
    return _info->status();
  }

  vector<string> Analysis::todos() const {
    if (!_info) return vector<string>();
    return _info->todos();
  }

  const vector<PdgIdPair>& Analysis::requiredBeams() const {
    return info().beams();
  }


  /// @todo Deprecate?
  Analysis& Analysis::setBeams(PdgId beam1, PdgId beam2) {
    assert(_info.get() != 0);
    _info->_beams.clear();
    _info->_beams += make_pair(beam1, beam2);
    return *this;
  }


  /// @todo Deprecate?
  bool Analysis::isCompatible(PdgId beam1, PdgId beam2) const {
    PdgIdPair beams(beam1, beam2);
    return isCompatible(beams);
  }


  /// @todo Deprecate?
  bool Analysis::isCompatible(const PdgIdPair& beams) const {
    foreach (const PdgIdPair& bp, requiredBeams()) {
      if (compatible(beams, bp)) return true;
    }
    return false;
    /// @todo Need to also check internal consistency of the analysis'
    /// beam requirements with those of the projections it uses.
  }


  Analysis& Analysis::setCrossSection(double xs) {
    _crossSection = xs;
    _gotCrossSection = true;
    return *this;
  }

  /// @todo Deprecate, eventually
  bool Analysis::needsCrossSection() const {
    return _needsCrossSection;
  }

  /// @todo Deprecate, eventually
  Analysis& Analysis::setNeedsCrossSection(bool needed) {
    _needsCrossSection = needed;
    return *this;
  }

  double Analysis::crossSection() const {
    if (!_gotCrossSection || _crossSection < 0) {
      string errMsg = "You did not set the cross section for the analysis " + name();
      throw Error(errMsg);
    }
    return _crossSection;
  }

  double Analysis::crossSectionPerEvent() const {
    const double sumW = sumOfWeights();
    assert(sumW > 0);
    return _crossSection / sumW;
  }


  AnalysisHandler& Analysis::handler() const {
    return *_analysishandler;
  }



  ////////////////////////////////////////////////////////////
  // Histogramming


  void Analysis::_cacheBinEdges() const {
    _cacheXAxisData();
    if (_histBinEdges.empty()) {
      getLog() << Log::TRACE << "Getting histo bin edges from AIDA for paper " << name() << endl;
      _histBinEdges = getBinEdges(_dpsData);
    }
  }


  void Analysis::_cacheXAxisData() const {
    if (_dpsData.empty()) {
      getLog() << Log::TRACE << "Getting DPS x-axis data from AIDA for paper " << name() << endl;
      _dpsData = getDPSXValsErrs(name());
    }
  }


  const BinEdges& Analysis::binEdges(const string& hname) const {
    _cacheBinEdges();
    getLog() << Log::TRACE << "Using histo bin edges for " << name() << ":" << hname << endl;
    const BinEdges& edges = _histBinEdges.find(hname)->second;
    if (getLog().isActive(Log::TRACE)) {
      stringstream edges_ss;
      foreach (const double be, edges) {
        edges_ss << " " << be;
      }
      getLog() << Log::TRACE << "Edges:" << edges_ss.str() << endl;
    }
    return edges;
  }


  const BinEdges& Analysis::binEdges(size_t datasetId, size_t xAxisId, size_t yAxisId) const {
    const string hname = makeAxisCode(datasetId, xAxisId, yAxisId);
    return binEdges(hname);
  }


  BinEdges Analysis::logBinEdges(size_t nbins, double lower, double upper) {
    assert(lower>0.0);
    assert(upper>lower);
    double loglower=log10(lower);
    double logupper=log10(upper);
    vector<double> binedges;
    double stepwidth=(logupper-loglower)/double(nbins);
    for (size_t i=0; i<=nbins; ++i) {
      binedges.push_back(pow(10.0, loglower+double(i)*stepwidth));
    }
    return binedges;
  }

  IHistogram1D* Analysis::bookHistogram1D(size_t datasetId, size_t xAxisId,
                                          size_t yAxisId, const string& title,
                                          const string& xtitle, const string& ytitle)
  {
    const string axisCode = makeAxisCode(datasetId, xAxisId, yAxisId);
    return bookHistogram1D(axisCode, title, xtitle, ytitle);
  }


  IHistogram1D* Analysis::bookHistogram1D(const string& hname, const string& title,
                                          const string& xtitle, const string& ytitle)
  {
    // Get the bin edges (only read the AIDA file once)
    const BinEdges edges = binEdges(hname);
    _makeHistoDir();
    const string path = histoPath(hname);
    IHistogram1D* hist = histogramFactory().createHistogram1D(path, title, edges);
    getLog() << Log::TRACE << "Made histogram " << hname <<  " for " << name() << endl;
    hist->setXTitle(xtitle);
    hist->setYTitle(ytitle);
    return hist;
  }


  IHistogram1D* Analysis::bookHistogram1D(const string& hname,
                                          size_t nbins, double lower, double upper,
                                          const string& title,
                                          const string& xtitle, const string& ytitle) {
    _makeHistoDir();
    const string path = histoPath(hname);
    IHistogram1D* hist = histogramFactory().createHistogram1D(path, title, nbins, lower, upper);
    getLog() << Log::TRACE << "Made histogram " << hname <<  " for " << name() << endl;
    hist->setXTitle(xtitle);
    hist->setYTitle(ytitle);
    return hist;
  }


  IHistogram1D* Analysis::bookHistogram1D(const string& hname,
                                          const vector<double>& binedges,
                                          const string& title,
                                          const string& xtitle, const string& ytitle) {
    _makeHistoDir();
    const string path = histoPath(hname);
    IHistogram1D* hist = histogramFactory().createHistogram1D(path, title, binedges);
    getLog() << Log::TRACE << "Made histogram " << hname <<  " for " << name() << endl;
    hist->setXTitle(xtitle);
    hist->setYTitle(ytitle);
    return hist;
  }


  /////////////////


  IProfile1D* Analysis::bookProfile1D(size_t datasetId, size_t xAxisId,
                                      size_t yAxisId, const string& title,
                                      const string& xtitle, const string& ytitle) {
    const string axisCode = makeAxisCode(datasetId, xAxisId, yAxisId);
    return bookProfile1D(axisCode, title, xtitle, ytitle);
  }


  IProfile1D* Analysis::bookProfile1D(const string& hname, const string& title,
                                      const string& xtitle, const string& ytitle)
  {
    // Get the bin edges (only read the AIDA file once)
    const BinEdges edges = binEdges(hname);
    _makeHistoDir();
    const string path = histoPath(hname);
    IProfile1D* prof = histogramFactory().createProfile1D(path, title, edges);
    getLog() << Log::TRACE << "Made profile histogram " << hname <<  " for " << name() << endl;
    prof->setXTitle(xtitle);
    prof->setYTitle(ytitle);
    return prof;
  }


  IProfile1D* Analysis::bookProfile1D(const string& hname,
                                      size_t nbins, double lower, double upper,
                                      const string& title,
                                      const string& xtitle, const string& ytitle) {
    _makeHistoDir();
    const string path = histoPath(hname);
    IProfile1D* prof = histogramFactory().createProfile1D(path, title, nbins, lower, upper);
    getLog() << Log::TRACE << "Made profile histogram " << hname <<  " for " << name() << endl;
    prof->setXTitle(xtitle);
    prof->setYTitle(ytitle);
    return prof;
  }


  IProfile1D* Analysis::bookProfile1D(const string& hname,
                                      const vector<double>& binedges,
                                      const string& title,
                                      const string& xtitle, const string& ytitle) {
    _makeHistoDir();
    const string path = histoPath(hname);
    IProfile1D* prof = histogramFactory().createProfile1D(path, title, binedges);
    getLog() << Log::TRACE << "Made profile histogram " << hname <<  " for " << name() << endl;
    prof->setXTitle(xtitle);
    prof->setYTitle(ytitle);
    return prof;
  }


  ///////////////////



  IDataPointSet* Analysis::bookDataPointSet(const string& hname, const string& title,
                                            const string& xtitle, const string& ytitle) {
    _makeHistoDir();
    const string path = histoPath(hname);
    IDataPointSet* dps = datapointsetFactory().create(path, title, 2);
    getLog() << Log::TRACE << "Made data point set " << hname <<  " for " << name() << endl;
    dps->setXTitle(xtitle);
    dps->setYTitle(ytitle);
    return dps;
  }


  IDataPointSet* Analysis::bookDataPointSet(const string& hname,
                                            size_t npts, double lower, double upper,
                                            const string& title,
                                            const string& xtitle, const string& ytitle) {
    IDataPointSet* dps = bookDataPointSet(hname, title, xtitle, ytitle);
    for (size_t pt = 0; pt < npts; ++pt) {
      const double binwidth = (upper-lower)/npts;
      const double bincentre = lower + (pt + 0.5) * binwidth;
      dps->addPoint();
      IMeasurement* meas = dps->point(pt)->coordinate(0);
      meas->setValue(bincentre);
      meas->setErrorPlus(binwidth/2.0);
      meas->setErrorMinus(binwidth/2.0);
    }
    return dps;
  }


  IDataPointSet* Analysis::bookDataPointSet(size_t datasetId, size_t xAxisId,
                                            size_t yAxisId, const string& title,
                                            const string& xtitle, const string& ytitle) {
    // Get the bin edges (only read the AIDA file once)
    _cacheXAxisData();
    // Build the axis code
    const string axisCode = makeAxisCode(datasetId, xAxisId, yAxisId);
    //const map<string, vector<DPSXPoint> > xpoints = getDPSXValsErrs(papername);
    getLog() << Log::TRACE << "Using DPS x-positions for " << name() << ":" << axisCode << endl;
    IDataPointSet* dps = bookDataPointSet(axisCode, title, xtitle, ytitle);
    const vector<DPSXPoint> xpts = _dpsData.find(axisCode)->second;
    for (size_t pt = 0; pt < xpts.size(); ++pt) {
      dps->addPoint();
      IMeasurement* meas = dps->point(pt)->coordinate(0);
      meas->setValue(xpts[pt].val);
      meas->setErrorPlus(xpts[pt].errplus);
      meas->setErrorMinus(xpts[pt].errminus);
    }
    getLog() << Log::TRACE << "Made DPS " << axisCode <<  " for " << name() << endl;
    return dps;
  }


  ////////////////////


  void Analysis::_makeHistoDir() {
    if (!_madeHistoDir) {
      if (! name().empty()) {
        // vector<string> dirs;
        // split(dirs, histoDir(), "/");
        // string pathpart;
        // foreach (const string& d, dirs) {
        //tree().mkdir();
        //}
        tree().mkdirs(histoDir());
      }
      _madeHistoDir = true;
    }
  }


  void Analysis::normalize(AIDA::IHistogram1D*& histo, double norm) {
    if (!histo) {
      getLog() << Log::ERROR << "Failed to normalise histo=NULL in analysis "
               << name() << " (norm=" << norm << ")" << endl;
      return;
    }
    const string hpath = tree().findPath(dynamic_cast<const AIDA::IManagedObject&>(*histo));
    getLog() << Log::TRACE << "Normalizing histo " << hpath << " to " << norm << endl;

    double oldintg = 0.0;
    int nBins = histo->axis().bins();
    for (int iBin = 0; iBin != nBins; ++iBin) {
      // Leaving out factor of binWidth because AIDA's "height" already includes a width factor.
      oldintg += histo->binHeight(iBin); // * histo->axis().binWidth(iBin);
    }
    if (oldintg == 0.0) {
      getLog() << Log::WARN << "Histo " << hpath << " has null integral during normalisation" << endl;
      return;
    }

    // Scale by the normalisation factor.
    scale(histo, norm/oldintg);
  }


  void Analysis::scale(AIDA::IHistogram1D*& histo, double scale) {
    if (!histo) {
      getLog() << Log::ERROR << "Failed to scale histo=NULL in analysis "
          << name() << " (scale=" << scale << ")" << endl;
      return;
    }
    const string hpath = tree().findPath(dynamic_cast<const AIDA::IManagedObject&>(*histo));
    getLog() << Log::TRACE << "Scaling histo " << hpath << endl;

    vector<double> x, y, ex, ey;
    for (size_t i = 0, N = histo->axis().bins(); i < N; ++i) {
      x.push_back(0.5 * (histo->axis().binLowerEdge(i) + histo->axis().binUpperEdge(i)));
      ex.push_back(histo->axis().binWidth(i)*0.5);

      // "Bin height" is a misnomer in the AIDA spec: width is neglected.
      // We'd like to do this: y.push_back(histo->binHeight(i) * scale);
      y.push_back(histo->binHeight(i)*scale/histo->axis().binWidth(i));

      // "Bin error" is a misnomer in the AIDA spec: width is neglected.
      // We'd like to do this: ey.push_back(histo->binError(i) * scale);
      ey.push_back(histo->binError(i)*scale/histo->axis().binWidth(i));
    }

    string title = histo->title();
    string xtitle = histo->xtitle();
    string ytitle = histo->ytitle();

    tree().mkdir("/tmpnormalize");
    tree().mv(hpath, "/tmpnormalize");

    AIDA::IDataPointSet* dps = datapointsetFactory().createXY(hpath, title, x, y, ex, ey);
    dps->setXTitle(xtitle);
    dps->setYTitle(ytitle);

    tree().rm(tree().findPath(dynamic_cast<AIDA::IManagedObject&>(*histo)));
    tree().rmdir("/tmpnormalize");

    // Set histo pointer to null - it can no longer be used.
    histo = 0;
  }


}
