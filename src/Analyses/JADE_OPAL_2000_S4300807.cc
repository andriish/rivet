// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/JADE_OPAL_2000_S4300807.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  // Constructor
  JADE_OPAL_2000_S4300807::JADE_OPAL_2000_S4300807(double sqrts, int nr_R_Jade,
						   int nr_R_Durham, int nr_y_Durham) :
    _sqrts(sqrts), _nr_R_Jade(nr_R_Jade),
    _nr_R_Durham(nr_R_Durham), _nr_y_Durham(nr_y_Durham)
  {
    setBeams(ELECTRON, POSITRON); 
    addProjection(Beam(), "Beams");
    const FinalState fs;
    addProjection(fs, "FS");
    addProjection(FastJets(fs, FastJets::JADE, 0.7), "JadeJets");
    addProjection(FastJets(fs, FastJets::DURHAM, 0.7), "DurhamJets");
  }


  void JADE_OPAL_2000_S4300807::analyze(const Event& e) {

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
      double y_23 = jadejet.clusterSeq()->exclusive_dmerge_max(2);
      double y_34 = jadejet.clusterSeq()->exclusive_dmerge_max(3);
      double y_45 = jadejet.clusterSeq()->exclusive_dmerge_max(4);
      double y_56 = jadejet.clusterSeq()->exclusive_dmerge_max(5);

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
      double y_23 = durjet.clusterSeq()->exclusive_dmerge_max(2);
      double y_34 = durjet.clusterSeq()->exclusive_dmerge_max(3);
      double y_45 = durjet.clusterSeq()->exclusive_dmerge_max(4);
      double y_56 = durjet.clusterSeq()->exclusive_dmerge_max(5);

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



  void JADE_OPAL_2000_S4300807::init() {
    stringstream ss;
    ss << _sqrts;
    _h_R_Jade[0]=bookDataPointSet(_nr_R_Jade, 1, 1, "Integrated 2-jet rate with Jade algorithm ("+ss.str()+" GeV)",
                                  "$y_{\\text{cut}}^\\text{Jade}$", "$R_2$");
    _h_R_Jade[1]=bookDataPointSet(_nr_R_Jade, 1, 2, "Integrated 3-jet rate with Jade algorithm ("+ss.str()+" GeV)",
                                  "$y_{\\text{cut}}^\\text{Jade}$", "$R_3$");
    _h_R_Jade[2]=bookDataPointSet(_nr_R_Jade, 1, 3, "Integrated 4-jet rate with Jade algorithm ("+ss.str()+" GeV)",
                                  "$y_{\\text{cut}}^\\text{Jade}$", "$R_4$");
    _h_R_Jade[3]=bookDataPointSet(_nr_R_Jade, 1, 4, "Integrated 5-jet rate with Jade algorithm ("+ss.str()+" GeV)",
                                  "$y_{\\text{cut}}^\\text{Jade}$", "$R_5$");
    _h_R_Jade[4]=bookDataPointSet(_nr_R_Jade, 1, 5, "Integrated $\\geq$6-jet rate with Jade algorithm ("+ss.str()+" GeV)",
                                  "$y_{\\text{cut}}^\\text{Jade}$", "$R_{\\geq 6}$");

    _h_R_Durham[0]=bookDataPointSet(_nr_R_Durham, 1, 1, "Integrated 2-jet rate with Durham algorithm ("+ss.str()+" GeV)",
                                    "$y_{\\text{cut}}^\\text{Durham}$", "$R_2$");
    _h_R_Durham[1]=bookDataPointSet(_nr_R_Durham, 1, 2, "Integrated 3-jet rate with Durham algorithm ("+ss.str()+" GeV)",
                                    "$y_{\\text{cut}}^\\text{Durham}$", "$R_3$");
    _h_R_Durham[2]=bookDataPointSet(_nr_R_Durham, 1, 3, "Integrated 4-jet rate with Durham algorithm ("+ss.str()+" GeV)",
                                    "$y_{\\text{cut}}^\\text{Durham}$", "$R_4$");
    _h_R_Durham[3]=bookDataPointSet(_nr_R_Durham, 1, 4, "Integrated 5-jet rate with Durham algorithm ("+ss.str()+" GeV)",
                                    "$y_{\\text{cut}}^\\text{Durham}$", "$R_5$");
    _h_R_Durham[4]=bookDataPointSet(_nr_R_Durham, 1, 5, "Integrated $\\geq$6-jet rate with Durham algorithm ("+ss.str()+" GeV)",
                                    "$y_{\\text{cut}}^\\text{Durham}$", "$R_{\\geq 6}$");

    _h_y_Durham[0]=bookHistogram1D(_nr_y_Durham, 1, 1, "Differential 2-jet rate with Durham algorithm ("+ss.str()+" GeV)",
                                   "$y_{23}^\\text{Durham}$", "$\\text{d}\\sigma/\\text{d}y_{23}$");
    _h_y_Durham[1]=bookHistogram1D(_nr_y_Durham, 1, 2, "Differential 3-jet rate with Durham algorithm ("+ss.str()+" GeV)",
                                   "$y_{34}^\\text{Durham}$", "$\\text{d}\\sigma/\\text{d}y_{34}$");
    _h_y_Durham[2]=bookHistogram1D(_nr_y_Durham, 1, 3, "Differential 4-jet rate with Durham algorithm ("+ss.str()+" GeV)",
                                   "$y_{45}^\\text{Durham}$", "$\\text{d}\\sigma/\\text{d}y_{45}$");
    _h_y_Durham[3]=bookHistogram1D(_nr_y_Durham, 1, 4, "Differential 5-jet rate with Durham algorithm ("+ss.str()+" GeV)",
                                   "$y_{56}^\\text{Durham}$", "$\\text{d}\\sigma/\\text{d}y_{56}$");
  }



  // Finalize
  void JADE_OPAL_2000_S4300807::finalize() {
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


}
