// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/AnalysisInfo.hh"
#include "Rivet/Tools/BeamConstraint.hh"

// #include "DummyConfig.hh"
// #ifdef HAVE_EXECINFO_H
// #include <execinfo.h>
// #endif

namespace Rivet {


  Analysis::Analysis(const string& name)
      : _analysishandler(NULL)
  {
    ProjectionApplier::_allowProjReg = false;
    _defaultname = name;

    unique_ptr<AnalysisInfo> ai = AnalysisInfo::make(name);
    assert(ai);
    _info = move(ai);
    assert(_info);
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
    /// @todo Cache in a member variable
    string _histoDir;
    if (_histoDir.empty()) {
      _histoDir = "/" + name();
      if (handler().runName().length() > 0) {
        _histoDir = "/" + handler().runName() + _histoDir;
      }
      replace_all(_histoDir, "//", "/"); //< iterates until none
    }
    return _histoDir;
  }


  const string Analysis::histoPath(const string& hname) const {
    const string path = histoDir() + "/" + hname;
    return path;
  }


  const string Analysis::histoPath(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) const {
    return histoDir() + "/" + mkAxisCode(datasetId, xAxisId, yAxisId);
  }


  const string Analysis::mkAxisCode(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) const {
    std::stringstream axisCode;
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


  Log& Analysis::getLog() const {
    string logname = "Rivet.Analysis." + name();
    return Log::getLog(logname);
  }


  ///////////////////////////////////////////


  size_t Analysis::numEvents() const {
    return handler().numEvents();
  }

  double Analysis::sumOfWeights() const {
    return handler().sumOfWeights();
  }


  ///////////////////////////////////////////


  bool Analysis::isCompatible(const ParticlePair& beams) const {
    return isCompatible(beams.first.pid(),  beams.second.pid(),
                        beams.first.energy(), beams.second.energy());
  }


  bool Analysis::isCompatible(PdgId beam1, PdgId beam2, double e1, double e2) const {
    PdgIdPair beams(beam1, beam2);
    pair<double,double> energies(e1, e2);
    return isCompatible(beams, energies);
  }


  bool Analysis::isCompatible(const PdgIdPair& beams, const pair<double,double>& energies) const {
    // First check the beam IDs
    bool beamIdsOk = false;
    for (const PdgIdPair& bp : requiredBeams()) {
      if (compatible(beams, bp)) {
        beamIdsOk =  true;
        break;
      }
    }
    if (!beamIdsOk) return false;

    // Next check that the energies are compatible (within 1% or 1 GeV, whichever is larger, for a bit of UI forgiveness)

    /// @todo Use some sort of standard ordering to improve comparisons, esp. when the two beams are different particles
    bool beamEnergiesOk = requiredEnergies().size() > 0 ? false : true;
    typedef pair<double,double> DoublePair;
    for (const DoublePair& ep : requiredEnergies()) {
      if ((fuzzyEquals(ep.first, energies.first, 0.01) && fuzzyEquals(ep.second, energies.second, 0.01)) ||
          (fuzzyEquals(ep.first, energies.second, 0.01) && fuzzyEquals(ep.second, energies.first, 0.01)) ||
          (abs(ep.first - energies.first) < 1*GeV && abs(ep.second - energies.second) < 1*GeV) ||
          (abs(ep.first - energies.second) < 1*GeV && abs(ep.second - energies.first) < 1*GeV)) {
        beamEnergiesOk =  true;
        break;
      }
    }
    return beamEnergiesOk;

    /// @todo Need to also check internal consistency of the analysis'
    /// beam requirements with those of the projections it uses.
  }


  ///////////////////////////////////////////

  double Analysis::crossSection() const {
    const YODA::Scatter1D::Points& ps = handler().crossSection()->points();
    if (ps.size() != 1) {
      string errMsg = "cross section missing for analysis " + name();
      throw Error(errMsg);
    }
    return ps[0].x();
  }

  double Analysis::crossSectionPerEvent() const {
    return crossSection()/sumOfWeights();
  }



  ////////////////////////////////////////////////////////////
  // Histogramming


  void Analysis::_cacheRefData() const {
    if (_refdata.empty()) {
      MSG_TRACE("Getting refdata cache for paper " << name());
      _refdata = getRefData(name());
    }
  }


  CounterPtr & Analysis::book(CounterPtr & ctr,
                      const string& cname) {
    const string path = histoPath(cname);
    ctr = CounterPtr(handler().weightNames(), Counter(path));
    addAnalysisObject(ctr);
    MSG_TRACE("Made counter " << cname << " for " << name());
    return ctr;
  }


  CounterPtr & Analysis::book(CounterPtr & ctr, unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) {
    const string axisCode = mkAxisCode(datasetId, xAxisId, yAxisId);
    return book(ctr, axisCode);
  }


  Histo1DPtr & Analysis::book(Histo1DPtr & histo, const string& hname, size_t nbins, double lower, double upper) {
    const string path = histoPath(hname);

    Histo1D hist = Histo1D(nbins, lower, upper, path);

    histo = Histo1DPtr(handler().weightNames(), hist);

    addAnalysisObject(histo);
    MSG_TRACE("Made histogram " << hname <<  " for " << name());
    return histo;
//    Histo1DPtr hist;
//    try { // try to bind to pre-existing
//      // AnalysisObjectPtr ao = getAnalysisObject(path);
//      // hist = dynamic_pointer_cast<Histo1D>(ao);
//      hist = getHisto1D(hname);
//      /// @todo Test that cast worked
//      /// @todo Also test that binning is as expected?
//      MSG_TRACE("Bound pre-existing histogram " << hname <<  " for " << name());
//    } catch (...) { // binding failed; make it from scratch
//      hist = make_shared<Histo1D>(nbins, lower, upper, histoPath(hname), title);
//      addAnalysisObject(hist);
//      MSG_TRACE("Made histogram " << hname <<  " for " << name());
//    }
  }

  Histo1DPtr & Analysis::book(Histo1DPtr & histo, const string& hname, const initializer_list<double>& binedges) {
  	return book(histo, hname, vector<double>{binedges});
  }

  Histo1DPtr & Analysis::book(Histo1DPtr & histo, const string& hname, const vector<double>& binedges) {
    const string path = histoPath(hname);
//    Histo1DPtr hist;
//    try { // try to bind to pre-existing
//      // AnalysisObjectPtr ao = getAnalysisObject(path);
//      // hist = dynamic_pointer_cast<Histo1D>(ao);
//      hist = getHisto1D(hname);
//      /// @todo Test that cast worked
//      /// @todo Also test that binning is as expected?
//      MSG_TRACE("Bound pre-existing histogram " << hname <<  " for " << name());
//    } catch (...) { // binding failed; make it from scratch
//      hist = make_shared<Histo1D>(binedges, histoPath(hname), title);
//      addAnalysisObject(hist);
//      MSG_TRACE("Made histogram " << hname <<  " for " << name());
//    }
    Histo1D hist = Histo1D(binedges, path);

    histo = Histo1DPtr(handler().weightNames(), hist);
    addAnalysisObject(histo);

    MSG_TRACE("Made histogram " << hname <<  " for " << name());
    return histo;
  }


  Histo1DPtr & Analysis::book(Histo1DPtr & histo, const string& hname, const Scatter2D& refscatter) {
    const string path = histoPath(hname);
//    Histo1DPtr hist;
//    try { // try to bind to pre-existing
//      // AnalysisObjectPtr ao = getAnalysisObject(path);
//     // hist = dynamic_pointer_cast<Histo1D>(ao);
//      hist = getHisto1D(hname);
//      /// @todo Test that cast worked
//      /// @todo Also test that binning is as expected?
//      MSG_TRACE("Bound pre-existing histogram " << hname <<  " for " << name());
//    } catch (...) { // binding failed; make it from scratch
//      hist = make_shared<Histo1D>(refscatter, histoPath(hname));
//      if (hist->hasAnnotation("IsRef")) hist->rmAnnotation("IsRef");
//      addAnalysisObject(hist);
//	MSG_TRACE("Made histogram " << hname <<  " for " << name());
//    }
    Histo1D hist = Histo1D(refscatter, path);

    histo = Histo1DPtr(handler().weightNames(), hist);
    addAnalysisObject(histo);

    MSG_TRACE("Made histogram " << hname <<  " for " << name());
    return histo;
  }

  Histo1DPtr & Analysis::book(Histo1DPtr & histo, const string& hname) {
    const Scatter2D& refdata = refData(hname);
    return book(histo, hname, refdata);
  }


  Histo1DPtr & Analysis::book(Histo1DPtr & histo, unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) {
    const string axisCode = mkAxisCode(datasetId, xAxisId, yAxisId);
    return book(histo, axisCode);
  }


  /// @todo Add booking methods which take a path, titles and *a reference Scatter from which to book*


  /////////////////


  Histo2DPtr & Analysis::book(Histo2DPtr & h2d,const string& hname,
                                   size_t nxbins, double xlower, double xupper,
                                   size_t nybins, double ylower, double yupper)
  {
    const string path = histoPath(hname);

    Histo2D hist(nxbins, xlower, xupper, nybins, ylower, yupper, path);

    h2d = Histo2DPtr(handler().weightNames(), hist);
    addAnalysisObject(h2d);

    MSG_TRACE("Made 2D histogram " << hname <<  " for " << name());
    return h2d;
  }


  Histo2DPtr & Analysis::book(Histo2DPtr & h2d,const string& hname,
                                   const initializer_list<double>& xbinedges,
                                   const initializer_list<double>& ybinedges)
  {
  	return book(h2d, hname, vector<double>{xbinedges}, vector<double>{ybinedges});
  }

  Histo2DPtr & Analysis::book(Histo2DPtr & h2d,const string& hname,
                                   const vector<double>& xbinedges,
                                   const vector<double>& ybinedges)
  {
    const string path = histoPath(hname);

    Histo2D hist(xbinedges, ybinedges, path);

    h2d = Histo2DPtr(handler().weightNames(), hist);
    addAnalysisObject(h2d);

    MSG_TRACE("Made 2D histogram " << hname <<  " for " << name());
    return h2d;
  }

  Histo2DPtr & Analysis::book(Histo2DPtr & histo, const string& hname,
                                   const Scatter3D& refscatter) {
    const string path = histoPath(hname);

    Histo2D hist = Histo2D(refscatter, path);

    histo = Histo2DPtr(handler().weightNames(), hist);
    addAnalysisObject(histo);

    MSG_TRACE("Made histogram " << hname <<  " for " << name());
    return histo;
  }

  Histo2DPtr & Analysis::book(Histo2DPtr & histo, const string& hname) {
    const Scatter3D& refdata = refData<Scatter3D>(hname);
    return book(histo, hname, refdata);
  }

  Histo2DPtr & Analysis::book(Histo2DPtr & histo, unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) {
    const string axisCode = mkAxisCode(datasetId, xAxisId, yAxisId);
    return book(histo, axisCode);
  }



  Profile1DPtr & Analysis::book(Profile1DPtr & p1d,const string& hname,
                                       size_t nbins, double lower, double upper) {
    const string path = histoPath(hname);

    Profile1D prof(nbins, lower, upper, path);

    p1d = Profile1DPtr(handler().weightNames(), prof);
    addAnalysisObject(p1d);

    MSG_TRACE("Made profile histogram " << hname <<  " for " << name());
    return p1d;
  }


  Profile1DPtr & Analysis::book(Profile1DPtr & p1d,const string& hname,
                                       const initializer_list<double>& binedges) {
  	return book(p1d, hname, vector<double>{binedges});
  }

  Profile1DPtr & Analysis::book(Profile1DPtr & p1d,const string& hname,
                                       const vector<double>& binedges) {
    const string path = histoPath(hname);
    Profile1D prof(binedges, path);

    p1d = Profile1DPtr(handler().weightNames(), prof);
    addAnalysisObject(p1d);

    MSG_TRACE("Made profile histogram " << hname <<  " for " << name());
    return p1d;
  }


  Profile1DPtr & Analysis::book(Profile1DPtr & p1d,const string& hname, const Scatter2D& refscatter) {
    const string path = histoPath(hname);
    Profile1D prof(refscatter, path);

    p1d = Profile1DPtr(handler().weightNames(), prof);
    addAnalysisObject(p1d);

    MSG_TRACE("Made profile histogram " << hname <<  " for " << name());
    return p1d;
    //if (prof.hasAnnotation("IsRef")) prof.rmAnnotation("IsRef");
  }


  Profile1DPtr & Analysis::book(Profile1DPtr & p1d,const string& hname) {
    const Scatter2D& refdata = refData(hname);
    book(p1d, hname, refdata);
    return p1d;
  }


  Profile1DPtr & Analysis::book(Profile1DPtr & p1d,unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) {
    const string axisCode = mkAxisCode(datasetId, xAxisId, yAxisId);
    return book(p1d, axisCode);
  }


  Profile2DPtr & Analysis::book(Profile2DPtr & p2d, const string& hname,
                                   size_t nxbins, double xlower, double xupper,
                                   size_t nybins, double ylower, double yupper)
  {
    const string path = histoPath(hname);
    Profile2D prof(nxbins, xlower, xupper, nybins, ylower, yupper, path);

    p2d = Profile2DPtr(handler().weightNames(), prof);
    addAnalysisObject(p2d);

    MSG_TRACE("Made 2D profile histogram " << hname <<  " for " << name());
    return p2d;
  }


  Profile2DPtr & Analysis::book(Profile2DPtr & p2d, const string& hname,
                                   const initializer_list<double>& xbinedges,
                                   const initializer_list<double>& ybinedges)
  {
  	return book(p2d, hname, vector<double>{xbinedges}, vector<double>{ybinedges});
  }

  Profile2DPtr & Analysis::book(Profile2DPtr & p2d, const string& hname,
                                   const vector<double>& xbinedges,
                                   const vector<double>& ybinedges)
  {
    const string path = histoPath(hname);
    Profile2D prof(xbinedges, ybinedges, path);

    p2d = Profile2DPtr(handler().weightNames(), prof);
    addAnalysisObject(p2d);

    MSG_TRACE("Made 2D profile histogram " << hname <<  " for " << name());
    return p2d;
  }


  Scatter2DPtr & Analysis::book(Scatter2DPtr & s2d, unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId,
                                       bool copy_pts) {
    const string axisCode = mkAxisCode(datasetId, xAxisId, yAxisId);
    return book(s2d, axisCode, copy_pts);
  }


  Scatter2DPtr & Analysis::book(Scatter2DPtr & s2d, const string& hname, bool copy_pts) {
    Scatter2D scat;
    const string path = histoPath(hname);
    if (copy_pts) {
      const Scatter2D& refdata = refData(hname);
      scat = Scatter2D(refdata, path);
      for (Point2D& p : scat.points()) p.setY(0, 0);
    } else {
      scat = Scatter2D(path);
    }

    s2d = Scatter2DPtr(handler().weightNames(), scat);
    addAnalysisObject(s2d);

    MSG_TRACE("Made scatter " << hname <<  " for " << name());
//    if (scat.hasAnnotation("IsRef")) scat.rmAnnotation("IsRef");
    return s2d;
  }


  Scatter2DPtr & Analysis::book(Scatter2DPtr & s2d, const string& hname,
                                       size_t npts, double lower, double upper) {
    // TODO: default branch has a read mechanism implemented, to start from an existing AO.
    // need to work out how to implement that for multiweights
    const string path = histoPath(hname);
    Scatter2D scat;
    const double binwidth = (upper-lower)/npts;
    for (size_t pt = 0; pt < npts; ++pt) {
      const double bincentre = lower + (pt + 0.5) * binwidth;
      scat.addPoint(bincentre, 0, 0.5*binwidth, 0);
    }

    s2d = Scatter2DPtr(handler().weightNames(), scat);
    addAnalysisObject(s2d);

    MSG_TRACE("Made scatter " << hname <<  " for " << name());
    return s2d;
  }


  Scatter2DPtr & Analysis::book(Scatter2DPtr & s2d, const string& hname, const vector<double>& binedges) {
    const string path = histoPath(hname);
    Scatter2D scat;
    for (size_t pt = 0; pt < binedges.size()-1; ++pt) {
      const double bincentre = (binedges[pt] + binedges[pt+1]) / 2.0;
      const double binwidth = binedges[pt+1] - binedges[pt];
      scat.addPoint(bincentre, 0, 0.5*binwidth, 0);
    }

    s2d = Scatter2DPtr(handler().weightNames(), scat);
    addAnalysisObject(s2d);

    MSG_TRACE("Made scatter " << hname <<  " for " << name());
    return s2d;
  }


  void Analysis::divide(CounterPtr c1, CounterPtr c2, Scatter1DPtr s) const {
    const string path = s->path();
    *s = *c1 / *c2;
    s->setPath(path);
  }

  void Analysis::divide(const Counter& c1, const Counter& c2, Scatter1DPtr s) const {
    const string path = s->path();
    *s = c1 / c2;
    s->setPath(path);
  }


  void Analysis::divide(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = *h1 / *h2;
    s->setPath(path);
  }

  void Analysis::divide(const Histo1D& h1, const Histo1D& h2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = h1 / h2;
    s->setPath(path);
  }


  void Analysis::divide(Profile1DPtr p1, Profile1DPtr p2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = *p1 / *p2;
    s->setPath(path);
  }

  void Analysis::divide(const Profile1D& p1, const Profile1D& p2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = p1 / p2;
    s->setPath(path);
  }


  void Analysis::divide(Histo2DPtr h1, Histo2DPtr h2, Scatter3DPtr s) const {
    const string path = s->path();
    *s = *h1 / *h2;
    s->setPath(path);
  }

  void Analysis::divide(const Histo2D& h1, const Histo2D& h2, Scatter3DPtr s) const {
    const string path = s->path();
    *s = h1 / h2;
    s->setPath(path);
  }


  void Analysis::divide(Profile2DPtr p1, Profile2DPtr p2, Scatter3DPtr s) const {
    const string path = s->path();
    *s = *p1 / *p2;
    s->setPath(path);
  }

  void Analysis::divide(const Profile2D& p1, const Profile2D& p2, Scatter3DPtr s) const {
    const string path = s->path();
    *s = p1 / p2;
    s->setPath(path);
  }


  /// @todo Counter and Histo2D efficiencies and asymms


  void Analysis::efficiency(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = YODA::efficiency(*h1, *h2);
    s->setPath(path);
  }

  void Analysis::efficiency(const Histo1D& h1, const Histo1D& h2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = YODA::efficiency(h1, h2);
    s->setPath(path);
  }


  void Analysis::asymm(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = YODA::asymm(*h1, *h2);
    s->setPath(path);
  }

  void Analysis::asymm(const Histo1D& h1, const Histo1D& h2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = YODA::asymm(h1, h2);
    s->setPath(path);
  }


  void Analysis::scale(CounterPtr cnt, Analysis::CounterAdapter factor) {
    if (!cnt) {
      MSG_WARNING("Failed to scale counter=NULL in analysis " << name() << " (scale=" << double(factor) << ")");
      return;
    }
    if (std::isnan(double(factor)) || std::isinf(double(factor))) {
      MSG_WARNING("Failed to scale counter=" << cnt->path() << " in analysis: " << name() << " (invalid scale factor = " << double(factor) << ")");
      factor = 0;
    }
    MSG_TRACE("Scaling counter " << cnt->path() << " by factor " << double(factor));
    try {
      cnt->scaleW(factor);
    } catch (YODA::Exception& we) {
      MSG_WARNING("Could not scale counter " << cnt->path());
      return;
    }
  }


  void Analysis::normalize(Histo1DPtr histo, Analysis::CounterAdapter norm, bool includeoverflows) {
    if (!histo) {
      MSG_WARNING("Failed to normalize histo=NULL in analysis " << name() << " (norm=" << double(norm) << ")");
      return;
    }
    MSG_TRACE("Normalizing histo " << histo->path() << " to " << double(norm));
    try {
      const double hint = histo->integral(includeoverflows);
      if (hint == 0)  MSG_WARNING("Skipping histo with null area " << histo->path());
      else            histo->normalize(norm, includeoverflows);
    } catch (YODA::Exception& we) {
      MSG_WARNING("Could not normalize histo " << histo->path());
      return;
    }
  }


  void Analysis::scale(Histo1DPtr histo, Analysis::CounterAdapter factor) {
    if (!histo) {
      MSG_WARNING("Failed to scale histo=NULL in analysis " << name() << " (scale=" << double(factor) << ")");
      return;
    }
    if (std::isnan(double(factor)) || std::isinf(double(factor))) {
      MSG_WARNING("Failed to scale histo=" << histo->path() << " in analysis: " << name() << " (invalid scale factor = " << double(factor) << ")");
      factor = 0;
    }
    MSG_TRACE("Scaling histo " << histo->path() << " by factor " << double(factor));
    try {
      histo->scaleW(factor);
    } catch (YODA::Exception& we) {
      MSG_WARNING("Could not scale histo " << histo->path());
      return;
    }
  }


  void Analysis::normalize(Histo2DPtr histo, Analysis::CounterAdapter norm, bool includeoverflows) {
    if (!histo) {
      MSG_ERROR("Failed to normalize histo=NULL in analysis " << name() << " (norm=" << double(norm) << ")");
      return;
    }
    MSG_TRACE("Normalizing histo " << histo->path() << " to " << double(norm));
    try {
      const double hint = histo->integral(includeoverflows);
      if (hint == 0)  MSG_WARNING("Skipping histo with null area " << histo->path());
      else            histo->normalize(norm, includeoverflows);
    } catch (YODA::Exception& we) {
      MSG_WARNING("Could not normalize histo " << histo->path());
      return;
    }
  }


  void Analysis::scale(Histo2DPtr histo, Analysis::CounterAdapter factor) {
    if (!histo) {
      MSG_ERROR("Failed to scale histo=NULL in analysis " << name() << " (scale=" << double(factor) << ")");
      return;
    }
    if (std::isnan(double(factor)) || std::isinf(double(factor))) {
      MSG_ERROR("Failed to scale histo=" << histo->path() << " in analysis: " << name() << " (invalid scale factor = " << double(factor) << ")");
      factor = 0;
    }
    MSG_TRACE("Scaling histo " << histo->path() << " by factor " << double(factor));
    try {
      histo->scaleW(factor);
    } catch (YODA::Exception& we) {
      MSG_WARNING("Could not scale histo " << histo->path());
      return;
    }
  }


  void Analysis::integrate(Histo1DPtr h, Scatter2DPtr s) const {
    // preserve the path info
    const string path = s->path();
    *s = toIntegralHisto(*h);
    s->setPath(path);
  }

  void Analysis::integrate(const Histo1D& h, Scatter2DPtr s) const {
    // preserve the path info
    const string path = s->path();
    *s = toIntegralHisto(h);
    s->setPath(path);
  }

}
  /// @todo 2D versions of integrate... defined how, exactly?!?


  //////////////////////////////////

namespace {
  void errormsg(std::string name) {
// #ifdef HAVE_BACKTRACE
//      void * buffer[4];
//      backtrace(buffer, 4);
//      backtrace_symbols_fd(buffer, 4 , 1);
// #endif
    std::cerr << name << ": Can't book objects outside of init().\n";
    assert(false);
  }
}

namespace Rivet {

  void Analysis::addAnalysisObject(const MultiweightAOPtr & ao) {
    if (handler().stage() == AnalysisHandler::Stage::INIT) {
      _analysisobjects.push_back(ao);
    }
    else {
      errormsg(name());
    }
  }


  void Analysis::removeAnalysisObject(const string& path) {
    for (auto it = _analysisobjects.begin();
         it != _analysisobjects.end(); ++it) {
      if ((*it)->path() == path) {
        _analysisobjects.erase(it);
        break;
      }
    }

  }

  void Analysis::removeAnalysisObject(const MultiweightAOPtr & ao) {
    for (auto it = _analysisobjects.begin();  it != _analysisobjects.end(); ++it) {
      if ((*it) == ao) {
        _analysisobjects.erase(it);
        break;
      }
    }
  }

}
