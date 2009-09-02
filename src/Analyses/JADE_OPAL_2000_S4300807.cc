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
    JADE_OPAL_2000_S4300807(const string& sqrtsstr, double sqrts, 
                            int nr_R_Jade, int nr_R_Durham, int nr_y_Durham)
      : Analysis("JADE_OPAL_2000_S4300807" + ("_" + sqrtsstr + "GEV")),
        _sqrts(sqrts), 
        _nr_R_Jade(nr_R_Jade),
        _nr_R_Durham(nr_R_Durham), 
        _nr_y_Durham(nr_y_Durham)
    {
      setBeams(ELECTRON, POSITRON); 
      addProjection(Beam(), "Beams");
      const FinalState fs;
      addProjection(fs, "FS");
      addProjection(FastJets(fs, FastJets::JADE, 0.7), "JadeJets");
      addProjection(FastJets(fs, FastJets::DURHAM, 0.7), "DurhamJets");
    }
    
    //@}


    /// @name Publication metadata
    //@{
    /// Collider on which the experiment ran.
    string collider() const {
      if (_sqrts < 90.0) {
        return "DESY PETRA";
      } else if (inRange(_sqrts, 90.0, 92.0)) {
        return "LEP Run I";
      } else {
        return "LEP Run 2";
      }
    }
    //@}

    
    /// @name Analysis methods
    //@{

    void init() {
      for (size_t i=0; i<5; ++i) {
        _h_R_Jade[i]=bookDataPointSet(_nr_R_Jade, 1, i+1);
        _h_R_Durham[i]=bookDataPointSet(_nr_R_Durham, 1, i+1);
        if (i<4)_h_y_Durham[i]=bookHistogram1D(_nr_y_Durham, 1, i+1);
      }
    }



    void analyze(const Event& e) {
      
      // Are we running with a compatible CMS energy?
      const double sbeams = applyProjection<Beam>(e, "Beams").sqrtS();
      if (fabs(sbeams - _sqrts)/GeV > 0.5) {
        getLog() << Log::ERROR 
                 << "CMS energy of events sqrt(s) = " << sbeams
                 <<" doesn't match analysis energy sqrt(s) = " << _sqrts << endl;
        /// @todo Really call exit()? I don't like the break of "command chain" that this implies
        exit(1);
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

    double _sqrts;
    int _nr_R_Jade, _nr_R_Durham, _nr_y_Durham;

  };



  //////////////////////////////////////////////////////////////



  class JADE_OPAL_2000_S4300807_35GEV : public JADE_OPAL_2000_S4300807 {
  public:
    JADE_OPAL_2000_S4300807_35GEV() : JADE_OPAL_2000_S4300807("35", 35.0, 7, 16, 24) {}
    string summary() const { return "Jet rates in e+e- at JADE [35 GeV]."; }
  };
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<JADE_OPAL_2000_S4300807_35GEV> plugin_JADE_OPAL_2000_S4300807_35GEV;


  class JADE_OPAL_2000_S4300807_44GEV : public JADE_OPAL_2000_S4300807 {
  public:
    JADE_OPAL_2000_S4300807_44GEV() : JADE_OPAL_2000_S4300807("44", 44.0, 8, 17, 25) {}
    string summary() const { return "Jet rates in e+e- at JADE [44 GeV]."; }
  };
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<JADE_OPAL_2000_S4300807_44GEV> plugin_JADE_OPAL_2000_S4300807_44GEV;


  class JADE_OPAL_2000_S4300807_91GEV : public JADE_OPAL_2000_S4300807 {
  public:
    JADE_OPAL_2000_S4300807_91GEV() : JADE_OPAL_2000_S4300807("91", 91.2, 9, 18, 26) {}
    string summary() const { return "Jet rates in e+e- at OPAL [91 GeV]."; }
  };
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<JADE_OPAL_2000_S4300807_91GEV> plugin_JADE_OPAL_2000_S4300807_91GEV;


  class JADE_OPAL_2000_S4300807_133GEV : public JADE_OPAL_2000_S4300807 {
  public:
    JADE_OPAL_2000_S4300807_133GEV() : JADE_OPAL_2000_S4300807("133", 133.0, 10, 19, 27) {}
    string summary() const { return "Jet rates in e+e- at OPAL [133 GeV]."; }
  };
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<JADE_OPAL_2000_S4300807_133GEV> plugin_JADE_OPAL_2000_S4300807_133GEV;


  class JADE_OPAL_2000_S4300807_161GEV : public JADE_OPAL_2000_S4300807 {
  public:
    JADE_OPAL_2000_S4300807_161GEV() : JADE_OPAL_2000_S4300807("161", 161.0, 11, 20, 28) {}
    string summary() const { return "Jet rates in e+e- at OPAL [161 GeV]."; }
  };
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<JADE_OPAL_2000_S4300807_161GEV> plugin_JADE_OPAL_2000_S4300807_161GEV;


  class JADE_OPAL_2000_S4300807_172GEV : public JADE_OPAL_2000_S4300807 {
  public:
    JADE_OPAL_2000_S4300807_172GEV() : JADE_OPAL_2000_S4300807("172", 172.0, 12, 21, 29) {}
    string summary() const { return "Jet rates in e+e- at OPAL [172 GeV]."; }
  };
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<JADE_OPAL_2000_S4300807_172GEV> plugin_JADE_OPAL_2000_S4300807_172GEV;


  class JADE_OPAL_2000_S4300807_183GEV : public JADE_OPAL_2000_S4300807 {
  public:
    JADE_OPAL_2000_S4300807_183GEV() : JADE_OPAL_2000_S4300807("183", 183.0, 13, 22, 30) {}
    string summary() const { return "Jet rates in e+e- at OPAL [183 GeV]."; }
  };
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<JADE_OPAL_2000_S4300807_183GEV> plugin_JADE_OPAL_2000_S4300807_183GEV;


  class JADE_OPAL_2000_S4300807_189GEV : public JADE_OPAL_2000_S4300807 {
  public:
    JADE_OPAL_2000_S4300807_189GEV() : JADE_OPAL_2000_S4300807("189", 189.0, 14, 23, 31) {}
    string summary() const { return "Jet rates in e+e- at OPAL [189 GeV]."; }
  };
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<JADE_OPAL_2000_S4300807_189GEV> plugin_JADE_OPAL_2000_S4300807_189GEV;


}
