// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "AIDA/IManagedObject.h"
using namespace AIDA;

namespace Rivet {


  Analysis::Analysis() 
    : _analysishandler(0),
      _madeHistoDir(false), _vetoedWeightSum(0)
  { 
    _gotCrossSection = false;
    setBeams(ANY, ANY);
    setNeedsCrossSection(false);
  }
  

  IAnalysisFactory& Analysis::analysisFactory() {
    return getHandler().analysisFactory();
  }


  ITree& Analysis::tree() {
    return getHandler().tree();
  }


  IHistogramFactory& Analysis::histogramFactory() {
    return getHandler().histogramFactory();
  }


  IDataPointSetFactory& Analysis::datapointsetFactory() {
    return getHandler().datapointsetFactory();
  }


  Log& Analysis::getLog() const {
    string logname = "Rivet.Analysis." + name();
    return Log::getLog(logname);
  }


  size_t Analysis::numEvents() const { return getHandler().numEvents(); }


  double Analysis::sumOfWeights() const { return getHandler().sumOfWeights() - _vetoedWeightSum; }


  void Analysis::_cacheBinEdges() {
    if (_histBinEdges.empty()) {
      getLog() << Log::TRACE << "Getting histo bin edges from AIDA for paper " << name() << endl;
      _histBinEdges = getBinEdges(name());
    }
  }


  string Analysis::_makeAxisCode(const size_t datasetId, const size_t xAxisId, const size_t yAxisId) const {
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


  IHistogram1D* Analysis::bookHistogram1D(const size_t datasetId, const size_t xAxisId, 
                                          const size_t yAxisId, const string& title) {
    const string axisCode = _makeAxisCode(datasetId, xAxisId, yAxisId);
    return bookHistogram1D(axisCode, title);
  }


  IHistogram1D* Analysis::bookHistogram1D(const string& hname, const string& title) {
    // Get the bin edges (only read the AIDA file once)
    _cacheBinEdges();
    getLog() << Log::TRACE << "Using histo bin edges for " << name() << ":" << hname << endl;
    const BinEdges edges = _histBinEdges.find(hname)->second;
    _makeHistoDir();
    const string path = histoDir() + "/" + hname;
    IHistogram1D* hist = histogramFactory().createHistogram1D(path, title, edges);
    getLog() << Log::TRACE << "Made histogram " << hname <<  " for " << name() << endl;
    return hist;
  }


  IHistogram1D* Analysis::bookHistogram1D(const string& hname, const string& title, 
                                          const size_t nbins, const double lower, const double upper) {
    _makeHistoDir();
    const string path = histoDir() + "/" + hname;
    return histogramFactory().createHistogram1D(path, title, nbins, lower, upper);
  }


  IHistogram1D* Analysis::bookHistogram1D(const string& hname, const string& title, 
                                          const vector<double>& binedges) {
    _makeHistoDir();
    const string path = histoDir() + "/" + hname;
    return histogramFactory().createHistogram1D(path, title, binedges);
  }



  IProfile1D* Analysis::bookProfile1D(const size_t datasetId, const size_t xAxisId, 
                                      const size_t yAxisId, const string& title) {
    const string axisCode = _makeAxisCode(datasetId, xAxisId, yAxisId);
    return bookProfile1D(axisCode, title);
  }


  IProfile1D* Analysis::bookProfile1D(const std::string& hname, const std::string& title) {
    Log& log = getLog();
    // Get the bin edges (only read the AIDA file once)
    _cacheBinEdges();
    log << Log::TRACE << "Using profile histo bin edges for " << name() << ":" << hname << endl;
    const BinEdges edges = _histBinEdges.find(hname)->second;
    if (log.isActive(Log::TRACE)) {
        stringstream edges_ss;
        foreach (const double be, edges) {
          edges_ss << " " << be;
        }
        log << Log::TRACE << "Edges:" << edges_ss.str() << endl;
    }
    _makeHistoDir();
    const string path = histoDir() + "/" + hname;
    IProfile1D* prof = histogramFactory().createProfile1D(path, title, edges);
    log << Log::TRACE << "Made profile histogram " << hname <<  " for " << name() << endl;
    return prof;
  }


  IProfile1D* Analysis::bookProfile1D(const string& hname, const string& title, 
                                      const size_t nbins, const double lower, const double upper) {
    _makeHistoDir();
    const string path = histoDir() + "/" + hname;
    return histogramFactory().createProfile1D(path, title, nbins, lower, upper);
  }


  IProfile1D* Analysis::bookProfile1D(const string& hname, const string& title, 
                                      const vector<double>& binedges) {
    _makeHistoDir();
    const string path = histoDir() + "/" + hname;
    return histogramFactory().createProfile1D(path, title, binedges);
  }


  /// @todo Convert to auto-book?
  IDataPointSet* Analysis::bookDataPointSet(const string& hname, const string& title) {
    _makeHistoDir();
    const string path = histoDir() + "/" + hname;
    return datapointsetFactory().create(path, title, 2);
  }


  IDataPointSet* Analysis::bookDataPointSet(const string& hname, const string& title, 
                                            const size_t npts, const double lower, const double upper) {
    IDataPointSet* dps = bookDataPointSet(hname, title);
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


  IDataPointSet* Analysis::bookDataPointSet(const size_t datasetId, const size_t xAxisId, 
                                            const size_t yAxisId, const string& title) {
    Log& log = getLog();
    // Get the bin edges (only read the AIDA file once)
    _cacheBinEdges();
    // Build the axis code
    const string axisCode = _makeAxisCode(datasetId, xAxisId, yAxisId);
    log << Log::TRACE << "Using DPS x-positions for " << name() << ":" << axisCode << endl;
    IDataPointSet* dps = bookDataPointSet(axisCode, title);
    const BinEdges edges = _histBinEdges.find(axisCode)->second;
    if (log.isActive(Log::TRACE)) {
      stringstream edges_ss;
      for (BinEdges::const_iterator be = edges.begin(); be != edges.end(); ++be) {
        edges_ss << " " << *be;
      }
      log << Log::TRACE << "Edges:" << edges_ss.str() << endl;
    }
    for (size_t pt = 0; pt < edges.size()-1; ++pt) {
      const double lower = edges[pt];
      const double upper = edges[pt+1];
      const double err = (upper-lower)/2.0;
      const double centre = (lower + upper)/2.0;
      dps->addPoint();
      IMeasurement* meas = dps->point(pt)->coordinate(0);
      meas->setValue(centre);
      meas->setErrorPlus(err);
      meas->setErrorMinus(err);
    }
    log << Log::TRACE << "Made DPS " << axisCode <<  " for " << name() << endl;
    return dps;
  }


  void Analysis::_makeHistoDir() {
    if (!_madeHistoDir) {
      if (! name().empty()) {
        tree().mkdir(histoDir());
      }
      _madeHistoDir = true;
    }
  }


  const Cuts Analysis::cuts() const {
    Cuts totalCuts = _cuts;
    foreach (ConstProjectionPtr p, getProjections()) {
      totalCuts.addCuts(p->cuts());
    }
    return totalCuts;
  }
  

  const bool Analysis::checkConsistency() const {
    // Check consistency of analysis beams with allowed beams of each contained projection.
    // set<ConstProjectionPtr> projections = getProjections();
    // for (set<ConstProjectionPtr>::const_iterator p = projections.begin(); p != projections.end(); ++p) {
    //   if (! compatible(getBeams(), (*p)->getBeamPairs()) ) {
    //     throw Error("Analysis " + name() + " beams are inconsistent with " 
    //                         + "allowed beams for projection " + (*p)->name());
    //   }
    // }
    // // Check the consistency of the accumulated cuts (throws if wrong).
    // getCuts().checkConsistency();
    return true;
  }

  
  void Analysis::normalize(AIDA::IHistogram1D*& histo, const double norm) {
    assert(histo);
    getLog() << Log::TRACE << "Normalizing histo " << histo->title() << " to " << norm << endl;
    
    double oldintg = 0.0;
    int nBins = histo->axis().bins();
    for (int iBin = 0; iBin != nBins; ++iBin) {
      // Leaving out factor of binWidth because AIDA's "height" already includes a width factor.
      oldintg += histo->binHeight(iBin); // * histo->axis().binWidth(iBin);
    }
    if (oldintg == 0.0) {
      getLog() << Log::WARN << "Histo " << histo->title() 
               << " has null integral during normalisation" << endl;
      return;
    }
  
    // Scale by the normalisation factor.
    scale(histo, norm/oldintg);
  }


  void Analysis::scale(AIDA::IHistogram1D*& histo, const double scale) {
    getLog() << Log::TRACE << "Scaling histo " << histo->title() << endl;
    
    std::vector<double> x, y, ex, ey;
    for (size_t i = 0, N = histo->axis().bins(); i < N; ++i) {
      x.push_back(0.5 * (histo->axis().binLowerEdge(i) + histo->axis().binUpperEdge(i)));
      ex.push_back(histo->axis().binWidth(i)*0.5);

      // "Bin height" is a misnomer in the AIDA spec: width is neglected.
      // We'd like to do this: y.push_back(histo->binHeight(i) * scale);
      y.push_back(histo->binHeight(i)*scale/histo->axis().binWidth(i));

      // "Bin error" is a misnomer in the AIDA spec: width is neglected.
      // We'd like to do this: ey.push_back(histo->binError(i) * scale);
      ey.push_back(histo->binError(i)*scale/(0.5*histo->axis().binWidth(i)));
    }
    
    std::string path =
      tree().findPath(dynamic_cast<AIDA::IManagedObject&>(*histo));
    std::string title = histo->title();
    tree().mkdir("/tmpnormalize");
    tree().mv(path, "/tmpnormalize");
    
    datapointsetFactory().createXY(path, title, x, y, ex, ey);
    
    tree().rm(tree().findPath(dynamic_cast<AIDA::IManagedObject&>(*histo)));
    tree().rmdir("/tmpnormalize");
    
    // Set histo pointer to null - it can no longer be used.
    histo = 0;
  }
  
  
}
