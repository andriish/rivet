// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief CDF Run II inclusive jet cross-section using the kT algorithm.
  /// @author James Monk
  class CDF_2007_S7057202 : public Analysis {
  public:

    /// Constructor
    CDF_2007_S7057202()
      : Analysis("CDF_2007_S7057202"),
        _minY(0.1), _maxY(0.7), _jetMinPT(54.0*GeV)
    {
      setBeams(PROTON, ANTIPROTON);
      //setSqrtS(1960*GeV);
      setNeedsCrossSection(true);
    }

 
    /// @name Analysis methods
    //@{

    /// Book histos and set counters for number of events passed in each one
    void init() {
      // Set up projections
      const FinalState fs;
      addProjection(FastJets(fs, FastJets::KT, 0.5), "JetsD05");
      addProjection(FastJets(fs, FastJets::KT, 0.7), "JetsD07");
      addProjection(FastJets(fs, FastJets::KT, 1.0), "JetsD10");

      // Book histos
      _histoD05 = bookHistogram1D(6, 1, 1);
      _histoD10 = bookHistogram1D(7, 1, 1);
      _binnedHistosD07.addHistogram(  0, 0.1, bookHistogram1D(1, 1, 1));
      _binnedHistosD07.addHistogram(0.1, 0.7, bookHistogram1D(2, 1, 1));
      _binnedHistosD07.addHistogram(0.7, 1.1, bookHistogram1D(3, 1, 1));
      _binnedHistosD07.addHistogram(1.1, 1.6, bookHistogram1D(4, 1, 1));
      _binnedHistosD07.addHistogram(1.6, 2.1, bookHistogram1D(5, 1, 1));
   
      size_t yind = 0;
      for (vector<AIDA::IHistogram1D*>::const_iterator histIt = _binnedHistosD07.getHistograms().begin();
           histIt != _binnedHistosD07.getHistograms().end(); ++histIt){
        _eventsPassed[*histIt] = 0.0;
        _yBinWidths[*histIt] = 2.0 * (_ybins[yind+1]-_ybins[yind]);
        ++yind;
      }
      _eventsPassed[_histoD05] = 0.0;
      _yBinWidths[_histoD05] = 2.0*(-_ybins[1]+_ybins[2]);
      _eventsPassed[_histoD10] = 0.0;
      _yBinWidths[_histoD10] = 2.0*(-_ybins[1]+_ybins[2]);
    }
 
 
    /// Do the analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
   
      const PseudoJets jetListD07 = applyProjection<FastJets>(event, "JetsD07").pseudoJets();
      set< IHistogram1D*> passed;
      /// @todo Use Jet interface rather than FastJet:PseudoJet
      for (PseudoJets::const_iterator jet = jetListD07.begin(); jet != jetListD07.end(); ++jet) {
        const double pt = jet->perp();
        if (pt > _jetMinPT) {
          AIDA::IHistogram1D* histo = _binnedHistosD07.fill(fabs(jet->rapidity()), pt, weight);
          if (histo != 0) {
            if (histo->coordToIndex(pt) != IAxis::OVERFLOW_BIN) {
              passed.insert(histo);
              _eventsPassed[histo] += weight;
            }
          }
        }
      }
   
      /// @todo Use Jet interface rather than FastJet:PseudoJet
      const PseudoJets jetListD05 = applyProjection<FastJets>(event, "JetsD05").pseudoJets();
      for (PseudoJets::const_iterator jet = jetListD05.begin(); jet != jetListD05.end(); ++jet) {
        const double pt = jet->perp();
        if (pt > _jetMinPT) {
          double rap = fabs(jet->rapidity());
          if (rap >= _minY && rap < _maxY){
            _histoD05->fill(pt, weight);
            if (_histoD05->coordToIndex(pt) != IAxis::OVERFLOW_BIN){
              passed.insert(_histoD05);
              _eventsPassed[_histoD05] += weight;
            }
          }
        }
      }
   
      /// @todo Use Jet interface rather than FastJet:PseudoJet
      const PseudoJets jetListD10 = applyProjection<FastJets>(event, "JetsD10").pseudoJets();
      for (PseudoJets::const_iterator jet = jetListD10.begin(); jet != jetListD10.end(); ++jet){
        const double pt = jet->perp();
        if (pt > _jetMinPT) {
          double rap = fabs(jet->rapidity());
          if (rap >= _minY && rap < _maxY){
            _histoD10->fill(pt, weight);
            if (_histoD10->coordToIndex(pt) != IAxis::OVERFLOW_BIN){
              passed.insert(_histoD10);
              _eventsPassed[_histoD10] += weight;
            }
          }
        }
      }
    }
 
 
    // Normalise histograms to cross-section
    void finalize() {
      const double xSecPerEvent = crossSection()/nanobarn / sumOfWeights();
      getLog() << Log::INFO << "Cross-section = " << crossSection()/nanobarn << " nb" << endl;
   
      for (map<IHistogram1D*,double>::iterator histIt = _eventsPassed.begin(),
             histJt = _yBinWidths.begin(); histIt != _eventsPassed.end(); ++histIt, ++histJt) {
        IHistogram1D* hist = histIt->first;
        const double xSec = xSecPerEvent * histIt->second / histJt->second;
        normalize(hist, xSec);
      }
    }
 
        //@}
 
  private:

    /// Rapidity range of histograms for R=0.05 and R=1 kt jets
    const double _minY, _maxY;
     
    /// Min jet \f$ p_T \f$ cut.
    /// @todo Make static const and UPPERCASE?
    const double _jetMinPT;
 
    /// Counter for the number of events analysed (actually the sum of weights, hence double).
    double _eventsTried;

    /// @name Histograms
    //@{
    /// The number of events in each histogram
    map<AIDA::IHistogram1D*, double> _eventsPassed;

    /// The y bin width of each histogram
    map<AIDA::IHistogram1D*, double> _yBinWidths;

    /// The y bin edge values
    static const double _ybins[6];

    /// Histograms in different eta regions
    BinnedHistogram<double> _binnedHistosD07;

    // Single histogram for the \f$R=0.5\f$ \f$k_\perp\f$ jets
    AIDA::IHistogram1D* _histoD05;

    // Single histogram for the \f$R=1.0\f$ \f$k_\perp\f$ jets
    AIDA::IHistogram1D* _histoD10;
    //@}

  };


  // Initialise static
  const double CDF_2007_S7057202::_ybins[] = { 0.0, 0.1, 0.7, 1.1, 1.6, 2.1 };


  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CDF_2007_S7057202> plugin_CDF_2007_S7057202;

}
