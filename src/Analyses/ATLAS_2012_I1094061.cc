// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class ATLAS_2012_I1094061 : public Analysis {

    ////////////////////////////////////////////////////////////////////////////
    /**
     * Little container to hold a pair of foreground and background histos and then
     * divide them at the end of the analysis
     */
    struct HistoPair{
      enum HistoType { FOREGROUND, BACKGROUND };

      HistoPair(): _analysis(0),
                   _h_foreground(0), _h_background(0), _d_final(0)
      {  }

      void init(int ds, int xaxis, int yaxis, ATLAS_2012_I1094061* analysis){
        _ds = ds;
        _xaxis = xaxis;
        _yaxis = yaxis;
        _analysis = analysis;
        ++HistoPair::_s_counter;
        const BinEdges& edges = _analysis->binEdges(_ds, _xaxis, _yaxis);
        string sCount = boost::lexical_cast<string>(HistoPair::_s_counter);
        _h_foreground = analysis->bookHistogram1D("tmpForeground" + sCount, edges);
        _h_background = analysis->bookHistogram1D("tmpBackground" + sCount, edges);
      }

      void fillForeground(double value, double weight){
        _h_foreground->fill(value, weight);
        _h_foreground->fill(-value, weight);
      }

      void fillBackground(double value, double weight){
        _h_background->fill(value, weight);
        _h_background->fill(-value, weight);
      }

      void fill(double value, double weight, HistoType type){

        switch(type){
          case FOREGROUND:
            fillForeground(value, weight);
            break;
          case BACKGROUND:
            fillBackground(value, weight);
            break;
        }
      }

      void finalize(double wgtSum, double bgWeight, double avNTracks){

        _h_foreground->scale(1. / wgtSum);
        _h_background->scale(1. / bgWeight);

        string histoPath = _analysis->histoPath(_ds, _xaxis, _yaxis);

        AIDA::IDataPointSet *final = _analysis->histogramFactory().divide(histoPath, *_h_foreground, *_h_background);

        for (int ii=0; ii!= final->size(); ++ii) {
          AIDA::IDataPoint* pt = final->point(ii);
          double y = pt->coordinate(1)->value();
          pt->coordinate(1)->setValue(y-(avNTracks - 1));
        }

        _analysis->histogramFactory().destroy(_h_foreground);
        _analysis->histogramFactory().destroy(_h_background);
      }

    private:

      int _ds, _xaxis, _yaxis;

      ATLAS_2012_I1094061 *_analysis;

      AIDA::IHistogram1D* _h_foreground;
      AIDA::IHistogram1D* _h_background;
      AIDA::IDataPointSet* _d_final;

      static short _s_counter;

    };


    ////////////////////////////////////////////////////////////////////////////

  public:

    ATLAS_2012_I1094061(): Analysis("ATLAS_2012_I1094061"),
    _minpT(100.*MeV), _etaMax(2.5), _nVersions(5), _version(0),
    _etaCut(2.), _phiCut(0.5*M_PI),
    _historyInclusive(_nVersions, ParticleVector()), _historyN20(_nVersions, ParticleVector()),
    _historyInclusiveWgts(_nVersions, 0.), _historyN20Wgts(_nVersions, 0.),
    _particleCountInclusive(0.), _particleCountN20(0.),
    _weightInclusive(0.), _weightN20(0.),
    _bgWeightInclusive(0.), _bgWeightN20(0.){

    }

    ////////////////////////////////////////////////////////////////////////////
    void init(){

      const ChargedFinalState cfs(-2.5, 2.5, _minpT);
      addProjection(cfs, "ChargedParticles");

      // Only do the multiplicity > 20 plots for 7 TeV collisions
      _doN20 = (fabs(sqrtS() - 7000.*GeV) < 0.1*GeV);

      int yaxis = (_doN20)? 2: 1;

      _hp_DEta_0_pi.init(1, 1, yaxis, this);
      _hp_DEta_0_pi2.init(2, 1, yaxis, this);
      _hp_DEta_pi2_pi.init(3, 1, yaxis, this);

      _hp_DPhi_0_2.init(4, 1, yaxis, this);
      _hp_DPhi_2_5.init(5, 1, yaxis, this);

      if(_doN20){

        yaxis = 3;

        _hp_N20_DEta_0_pi.init(1, 1, yaxis, this);
        _hp_N20_DEta_0_pi2.init(2, 1, yaxis, this);
        _hp_N20_DEta_pi2_pi.init(3, 1, yaxis, this);

        _hp_N20_DPhi_0_2.init(4, 1, yaxis, this);
        _hp_N20_DPhi_2_5.init(5, 1, yaxis, this);

      }
      return;
    }

    ////////////////////////////////////////////////////////////////////////////
    void analyze(const Event &evt){

      const ChargedFinalState &cfsProj = applyProjection<ChargedFinalState>(evt, "ChargedParticles");

      ParticleVector chargedParticles = cfsProj.particles();

      if(chargedParticles.size() < 2) vetoEvent;

      bool hasN20 = (_doN20 && chargedParticles.size() >= 20);

      double dMultiplicity = (double)chargedParticles.size();

      double multiplicityWeightIncr = dMultiplicity * evt.weight();

      _weightInclusive += evt.weight();
      _particleCountInclusive += multiplicityWeightIncr;

      if(hasN20){
        _weightN20 += evt.weight();
        _particleCountN20 += multiplicityWeightIncr;
      }

      double fgWeight = 2.*evt.weight() / dMultiplicity;

      for(ParticleVector::const_iterator p1 = chargedParticles.begin();
          p1 != chargedParticles.end(); ++p1){

        ParticleVector::const_iterator p2 = p1;
        ++p2;

        // fill the foreground distributions
        while(p2 != chargedParticles.end()){
          fillHistosInclusive(*p1, *p2, fgWeight, HistoPair::FOREGROUND);
          if(hasN20) fillHistosN20(*p1, *p2, fgWeight, HistoPair::FOREGROUND);
          ++p2;
        }// end filling the foreground distributions

        // loop over the history of particles from previous events and fill the background
        // by correlating those particles with the current event

        for(size_t version = 0; version != _nVersions; ++version){

          const ParticleVector &bgParticles = _historyInclusive[version];
          double bgWeight = evt.weight() * _historyInclusiveWgts[version];

          for(ParticleVector::const_iterator p2 = bgParticles.begin();
              p2 != bgParticles.end(); ++p2){
            fillHistosInclusive(*p1, *p2, bgWeight, HistoPair::BACKGROUND);
            _bgWeightInclusive += bgWeight;
          }

          if(!hasN20) continue;

          const ParticleVector &bgParticlesN20 = _historyN20[version];
          bgWeight = evt.weight() * _historyN20Wgts[version];

          for(ParticleVector::const_iterator p2 = bgParticlesN20.begin();
              p2 != bgParticlesN20.end(); ++p2){
            fillHistosN20(*p1, *p2, bgWeight, HistoPair::BACKGROUND);
            _bgWeightN20 += bgWeight;
          }

        }//end loop over particle history for background fill

      }// end particle loop

      // Overwrite the history for the version count number
      _historyInclusive[_version] = chargedParticles;
      _historyInclusiveWgts[_version] = evt.weight();

      if(hasN20){
        _historyN20[_version] = chargedParticles;
        _historyN20Wgts[_version] = evt.weight();
      }

      ++_version;
      if(_version == _nVersions) _version = 0;

      return;
    }

    ////////////////////////////////////////////////////////////////////////////
    void finalize(){

      double avMultiplicity = _particleCountInclusive / _weightInclusive;

      _hp_DEta_0_pi.finalize(_weightInclusive,  _bgWeightInclusive, avMultiplicity);
      _hp_DEta_0_pi2.finalize(_weightInclusive, _bgWeightInclusive, avMultiplicity);
      _hp_DEta_pi2_pi.finalize(_weightInclusive,_bgWeightInclusive, avMultiplicity);

      _hp_DPhi_0_2.finalize(_weightInclusive, _bgWeightInclusive, avMultiplicity);
      _hp_DPhi_2_5.finalize(_weightInclusive, _bgWeightInclusive, avMultiplicity);

      if(_doN20){
        avMultiplicity = _particleCountN20 / _weightN20;
        _hp_N20_DEta_0_pi.finalize(_weightN20,   _bgWeightN20, avMultiplicity);
        _hp_N20_DEta_0_pi2.finalize(_weightN20,  _bgWeightN20, avMultiplicity);
        _hp_N20_DEta_pi2_pi.finalize(_weightN20, _bgWeightN20, avMultiplicity);

        _hp_N20_DPhi_0_2.finalize(_weightN20, _bgWeightN20, avMultiplicity);
        _hp_N20_DPhi_2_5.finalize(_weightN20, _bgWeightN20, avMultiplicity);
      }

      return;
    }

    ////////////////////////////////////////////////////////////////////////////

    void fillHistos(const Particle &p1, const Particle &p2, double weight,
                    HistoPair::HistoType type, bool inclusive){

      double dEta = fabs(p1.eta() - p2.eta());
      double dPhi = mapAngle0ToPi(p1.phi() - p2.phi());
      double dPhiShift = TWOPI - dPhi;

      HistoPair &dEta_0_pi   = (inclusive)? _hp_DEta_0_pi   :_hp_N20_DEta_0_pi;
      HistoPair &dPhi_0_2    = (inclusive)? _hp_DPhi_0_2    :_hp_N20_DPhi_0_2;
      HistoPair &dPhi_2_5    = (inclusive)? _hp_DPhi_2_5    :_hp_N20_DPhi_2_5;
      HistoPair &dEta_0_pi2  = (inclusive)? _hp_DEta_0_pi2  :_hp_N20_DEta_0_pi2;
      HistoPair &dEta_pi2_pi = (inclusive)? _hp_DEta_pi2_pi :_hp_N20_DEta_pi2_pi;

      dEta_0_pi.fill(dEta, weight, type);

      if(dEta < _etaCut){
        dPhi_0_2.fill(dPhi, weight, type);
        dPhi_0_2.fill(dPhiShift, weight, type);
      }else{
        dPhi_2_5.fill(dPhi, weight, type);
        dPhi_2_5.fill(dPhiShift, weight, type);
      }

      if(dPhi < _phiCut){
        dEta_0_pi2.fill(dEta, weight, type);
      }else{
        dEta_pi2_pi.fill(dEta, weight, type);
      }

      return;

    }
    ////////////////////////////////////////////////////////////////////////////

    void fillHistosInclusive(const Particle &p1, const Particle &p2, double weight,
                             HistoPair::HistoType type){

      fillHistos(p1, p2, weight, type, true);
      return;
    }

    void fillHistosN20(const Particle &p1, const Particle &p2, double weight,
                       HistoPair::HistoType type){

      fillHistos(p1, p2, weight, type, false);
      return;
    }
    ////////////////////////////////////////////////////////////////////////////


    double _minpT;
    double _etaMax;

    size_t _nVersions;
    size_t _version;

    double _etaCut;
    double _phiCut;

    /// The "history" vectors contain the history of particles from _nVersions previous events
    /// These are used to construct the background correlation.
    vector<ParticleVector> _historyInclusive;
    vector<ParticleVector> _historyN20;

    vector<double>         _historyInclusiveWgts;
    vector<double>         _historyN20Wgts;

    double _particleCountInclusive;
    double _particleCountN20;
    double _weightInclusive;
    double _weightN20;

    double _bgWeightInclusive;
    double _bgWeightN20;

    bool _doN20;

    HistoPair _hp_DEta_0_pi;
    HistoPair _hp_DEta_0_pi2;
    HistoPair _hp_DEta_pi2_pi;

    HistoPair _hp_DPhi_0_2;
    HistoPair _hp_DPhi_2_5;

    HistoPair _hp_N20_DEta_0_pi;
    HistoPair _hp_N20_DEta_0_pi2;
    HistoPair _hp_N20_DEta_pi2_pi;

    HistoPair _hp_N20_DPhi_0_2;
    HistoPair _hp_N20_DPhi_2_5;

  };


  short ATLAS_2012_I1094061::HistoPair::_s_counter = 0;


  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1094061);

}
