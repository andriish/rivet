// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /**
   * @brief Jet rates in e+e- at OPAL and JADE
   * @author Frank Siegert
   *
   * @par Run conditions
   *
   * @arg LEP1 beam energy: \f$ \sqrt{s} = \$f 91.2 GeV
   * @arg Run with generic QCD events.
   */
  class JADE_OPAL_2000_S4300807 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    JADE_OPAL_2000_S4300807() : Analysis("JADE_OPAL_2000_S4300807"),
        _initialised(false)
    {
      setBeams(ELECTRON, POSITRON); 
    }
    
    //@}

    
    /// @name Analysis methods
    //@{

    void init() {
      addProjection(Beam(), "Beams");
      const FinalState fs;
      addProjection(fs, "FS");
      addProjection(FastJets(fs, FastJets::JADE, 0.7), "JadeJets");
      addProjection(FastJets(fs, FastJets::DURHAM, 0.7), "DurhamJets");
    }



    void analyze(const Event& e) {
      
      // Which CMS energy are we running at?
      if (!_initialised) {
        const double sqrts = applyProjection<Beam>(e, "Beams").sqrtS()/GeV;
        int offset(0);
        switch (int(sqrts+0.5)) {
        case 35: offset = 7; break;
        case 44: offset = 8; break;
        case 91: offset = 9; break;
        case 133: offset = 10; break;
        case 161: offset = 11; break;
        case 172: offset = 12; break;
        case 183: offset = 13; break;
        case 189: offset = 14; break;
        default:
          getLog() << Log::ERROR 
              << "CMS energy of events sqrt(s) = " << sqrts
              <<" doesn't match any available analysis energy." << endl;
          /// @todo Really call exit()? I don't like the break of "command chain" that this implies
          exit(1);
        }
        for (size_t i=0; i<5; ++i) {
          _h_R_Jade[i]=bookDataPointSet(offset, 1, i+1);
          _h_R_Durham[i]=bookDataPointSet(offset+9, 1, i+1);
          if (i<4)_h_y_Durham[i]=bookHistogram1D(offset+17, 1, i+1);
        }
        _initialised = true;
      }
        
      
      // Jets
      getLog() << Log::DEBUG << "Using FastJet JADE patch to make diff jet rate plots:" << endl;
      const double weight = e.weight();
      
      const FastJets& jadejet = applyProjection<FastJets>(e, "JadeJets");
      if (jadejet.clusterSeq()) {
        double y_23 = jadejet.clusterSeq()->exclusive_ymerge(2);
        double y_34 = jadejet.clusterSeq()->exclusive_ymerge(3);
        double y_45 = jadejet.clusterSeq()->exclusive_ymerge(4);
        double y_56 = jadejet.clusterSeq()->exclusive_ymerge(5);
        
        for (int i = 0; i < _h_R_Jade[0]->size(); ++i) {
          IDataPoint* dp = _h_R_Jade[0]->point(i);
          if (y_23 < dp->coordinate(0)->value()) {
            dp->coordinate(1)->setValue(dp->coordinate(1)->value() + weight);
          }
        }
        for (int i = 0; i < _h_R_Jade[1]->size(); ++i) {
          IDataPoint* dp = _h_R_Jade[1]->point(i);
          double ycut = dp->coordinate(0)->value();
          if (y_34 < ycut && y_23 > ycut) {
            dp->coordinate(1)->setValue(dp->coordinate(1)->value() + weight);
          }
        }
        for (int i = 0; i < _h_R_Jade[2]->size(); ++i) {
          IDataPoint* dp = _h_R_Jade[2]->point(i);
          double ycut = dp->coordinate(0)->value();
          if (y_45 < ycut && y_34 > ycut) {
            dp->coordinate(1)->setValue(dp->coordinate(1)->value() + weight);
          }
        }
        for (int i = 0; i < _h_R_Jade[3]->size(); ++i) {
          IDataPoint* dp = _h_R_Jade[3]->point(i);
          double ycut = dp->coordinate(0)->value();
          if (y_56 < ycut && y_45 > ycut) {
            dp->coordinate(1)->setValue(dp->coordinate(1)->value() + weight);
          }
        }
        for (int i = 0; i < _h_R_Jade[4]->size(); ++i) {
          IDataPoint* dp = _h_R_Jade[4]->point(i);
          double ycut = dp->coordinate(0)->value();
          if (y_56 > ycut) {
            dp->coordinate(1)->setValue(dp->coordinate(1)->value() + weight);
          }
        }
      }
      
      const FastJets& durjet = applyProjection<FastJets>(e, "DurhamJets");
      if (durjet.clusterSeq()) {
        double y_23 = durjet.clusterSeq()->exclusive_ymerge(2);
        double y_34 = durjet.clusterSeq()->exclusive_ymerge(3);
        double y_45 = durjet.clusterSeq()->exclusive_ymerge(4);
        double y_56 = durjet.clusterSeq()->exclusive_ymerge(5);
        
        _h_y_Durham[0]->fill(y_23, weight);
        _h_y_Durham[1]->fill(y_34, weight);
        _h_y_Durham[2]->fill(y_45, weight);
        _h_y_Durham[3]->fill(y_56, weight);
        
        for (int i = 0; i < _h_R_Durham[0]->size(); ++i) {
          IDataPoint* dp = _h_R_Durham[0]->point(i);
          if (y_23 < dp->coordinate(0)->value()) {
            dp->coordinate(1)->setValue(dp->coordinate(1)->value() + weight);
          }
        }
        for (int i = 0; i < _h_R_Durham[1]->size(); ++i) {
          IDataPoint* dp = _h_R_Durham[1]->point(i);
          double ycut = dp->coordinate(0)->value();
          if (y_34 < ycut && y_23 > ycut) {
            dp->coordinate(1)->setValue(dp->coordinate(1)->value() + weight);
          }
        }
        for (int i = 0; i < _h_R_Durham[2]->size(); ++i) {
          IDataPoint* dp = _h_R_Durham[2]->point(i);
          double ycut = dp->coordinate(0)->value();
          if (y_45 < ycut && y_34 > ycut) {
            dp->coordinate(1)->setValue(dp->coordinate(1)->value() + weight);
          }
        }
        for (int i = 0; i < _h_R_Durham[3]->size(); ++i) {
          IDataPoint* dp = _h_R_Durham[3]->point(i);
          double ycut = dp->coordinate(0)->value();
          if (y_56 < ycut && y_45 > ycut) {
            dp->coordinate(1)->setValue(dp->coordinate(1)->value() + weight);
          }
        }
        for (int i = 0; i < _h_R_Durham[4]->size(); ++i) {
          IDataPoint* dp = _h_R_Durham[4]->point(i);
          double ycut = dp->coordinate(0)->value();
          if (y_56 > ycut) {
            dp->coordinate(1)->setValue(dp->coordinate(1)->value() + weight);
          }
        }
      }
    }



    /// Finalize
    void finalize() {
      for (size_t n = 0; n < 4; ++n) {
        scale(_h_y_Durham[n], 1.0/sumOfWeights());
      }
      
      for (size_t n = 0; n < 5; ++n) {
        /// scale integrated jet rates to 100%
        for (int i = 0; i < _h_R_Jade[n]->size(); ++i) {
          IDataPoint* dp = _h_R_Jade[n]->point(i);
          dp->coordinate(1)->setValue(dp->coordinate(1)->value()*100.0/sumOfWeights());
        }
        for (int i = 0; i < _h_R_Durham[n]->size(); ++i) {
          IDataPoint* dp = _h_R_Durham[n]->point(i);
          dp->coordinate(1)->setValue(dp->coordinate(1)->value()*100.0/sumOfWeights());
        }
      }
    }
    
    //@}
    
    
  private:

    /// @name Histograms
    //@{
    AIDA::IDataPointSet *_h_R_Jade[5];
    AIDA::IDataPointSet *_h_R_Durham[5];
    AIDA::IHistogram1D *_h_y_Durham[4];
    //@}

    bool _initialised;

  };



  //////////////////////////////////////////////////////////////



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<JADE_OPAL_2000_S4300807> plugin_JADE_OPAL_2000_S4300807;


}
