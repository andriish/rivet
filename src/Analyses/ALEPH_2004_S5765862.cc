// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"

namespace Rivet {


  class ALEPH_2004_S5765862 : public Analysis {
  public:

    ALEPH_2004_S5765862()
      : Analysis("ALEPH_2004_S5765862") , _initialised(false)
    {
      setBeams(ELECTRON, POSITRON);
    }


  public:

    void init() {
      _initialised=true;

      const FinalState fs;
      addProjection(fs, "FS");
      addProjection(FastJets(fs, FastJets::DURHAM, 0.7), "DurhamJets");

      const Thrust thrust(fs);
      addProjection(thrust, "Thrust");
      addProjection(Sphericity(fs), "Sphericity");
      addProjection(ParisiTensor(fs), "Parisi");
      addProjection(Hemispheres(thrust), "Hemispheres");
      
      // Histos
      int offset = 0;
      switch (int(sqrtS()/GeV + 0.5)) {
      case 91: offset = 0; break;
      case 133: offset = 1; break;
      case 161: offset = 2; break;
      case 172: offset = 3; break;
      case 183: offset = 4; break;
      case 189: offset = 5; break;
      case 200: offset = 6; break;
      case 206: offset = 7; break;
      default:
        getLog() << Log::WARNING
                 << "CMS energy of events sqrt(s) = " << sqrtS()/GeV
                 <<" doesn't match any available analysis energy." << endl;
        _initialised=false;
        return;
      }
      
      // event shapes
      _h_thrust = bookHistogram1D(offset+54, 1, 1);
      _h_heavyjetmass = bookHistogram1D(offset+62, 1, 1);
      _h_totaljetbroadening = bookHistogram1D(offset+70, 1, 1);
      _h_widejetbroadening = bookHistogram1D(offset+78, 1, 1);
      _h_cparameter = bookHistogram1D(offset+86, 1, 1);
      _h_thrustmajor = bookHistogram1D(offset+94, 1, 1);
      _h_thrustminor = bookHistogram1D(offset+102, 1, 1);
      _h_jetmassdifference = bookHistogram1D(offset+110, 1, 1);
      _h_aplanarity = bookHistogram1D(offset+118, 1, 1);
      // planarity is missing the 91 gev data, so left out
      _h_oblateness = bookHistogram1D(offset+133, 1, 1);
      _h_sphericity = bookHistogram1D(offset+141, 1, 1);
      
      
      // Durham n->m jet resolutions
      _h_y_Durham[0] = bookHistogram1D(offset+149, 1, 1);   // y12 d149 ... d156
      _h_y_Durham[1] = bookHistogram1D(offset+157, 1, 1);   // y23 d157 ... d164
      _h_y_Durham[2] = bookHistogram1D(offset+165, 1, 1);   // y34 d165 ... d172
      if (offset<6) { // there is no y45 and y56 for 200 gev
        _h_y_Durham[3] = bookHistogram1D(offset+173, 1, 1); // y45 d173 ... d179
        _h_y_Durham[4] = bookHistogram1D(offset+180, 1, 1); // y56 d180 ... d186
      }
      else if (offset==6) {
        _h_y_Durham[3] = NULL;
        _h_y_Durham[4] = NULL;
      }
      else if (offset==7) {
        _h_y_Durham[3] = bookHistogram1D(179, 1, 1);
        _h_y_Durham[4] = bookHistogram1D(186, 1, 1);
      }
      
      // Durham n-jet fractions      
      _h_R_Durham[0] = bookDataPointSet(offset+187, 1, 1); // R1 d187 ... d194
      _h_R_Durham[1] = bookDataPointSet(offset+195, 1, 1); // R2 d195 ... d202
      _h_R_Durham[2] = bookDataPointSet(offset+203, 1, 1); // R3 d203 ... d210
      _h_R_Durham[3] = bookDataPointSet(offset+211, 1, 1); // R4 d211 ... d218
      _h_R_Durham[4] = bookDataPointSet(offset+219, 1, 1); // R5 d219 ... d226
      _h_R_Durham[5] = bookDataPointSet(offset+227, 1, 1); // R>=6 d227 ... d234
    }


    void analyze(const Event& e) {
      if (!_initialised) return;
      const double weight = e.weight();

      // event shapes
      const Thrust& thrust = applyProjection<Thrust>(e, "Thrust");
      double thr = (fuzzyEquals(sqrtS(),91.2*GeV,0.5)?thrust.thrust():1.0-thrust.thrust());
      _h_thrust->fill(thr,weight);
      _h_thrustmajor->fill(thrust.thrustMajor(),weight);
      _h_thrustminor->fill(log(thrust.thrustMinor()),weight);
      _h_oblateness->fill(thrust.oblateness(),weight);
      
      const Hemispheres& hemi = applyProjection<Hemispheres>(e, "Hemispheres");
      _h_heavyjetmass->fill(hemi.scaledM2high(),weight);
      _h_jetmassdifference->fill(hemi.scaledM2diff(),weight);
      _h_totaljetbroadening->fill(hemi.Bsum(),weight);
      _h_widejetbroadening->fill(hemi.Bmax(),weight);
      
      const ParisiTensor& parisi = applyProjection<ParisiTensor>(e, "Parisi");
      _h_cparameter->fill(parisi.C(),weight);

      const Sphericity& sphericity = applyProjection<Sphericity>(e, "Sphericity");
      _h_aplanarity->fill(sphericity.aplanarity(),weight);
      _h_sphericity->fill(sphericity.sphericity(),weight);

      // jet rates
      const FastJets& durjet = applyProjection<FastJets>(e, "DurhamJets");
      if (durjet.clusterSeq()) {
        /// @todo Put this in an index loop?
        double y_12 = log(durjet.clusterSeq()->exclusive_ymerge(1));
        double y_23 = log(durjet.clusterSeq()->exclusive_ymerge(2));
        double y_34 = log(durjet.clusterSeq()->exclusive_ymerge(3));
        double y_45 = log(durjet.clusterSeq()->exclusive_ymerge(4));
        double y_56 = log(durjet.clusterSeq()->exclusive_ymerge(5));
     
        _h_y_Durham[0]->fill(-y_12, weight);
        _h_y_Durham[1]->fill(-y_23, weight);
        _h_y_Durham[2]->fill(-y_34, weight);
        _h_y_Durham[3]->fill(-y_45, weight);
        _h_y_Durham[4]->fill(-y_56, weight);
     
        /// @todo Do this more elegant, maybe even in finalize
        for (int i = 0; i < _h_R_Durham[0]->size(); ++i) {
          IDataPoint* dp = _h_R_Durham[0]->point(i);
          if (y_12 < dp->coordinate(0)->value()) {
            dp->coordinate(1)->setValue(dp->coordinate(1)->value() + weight);
          }
        }
        for (int i = 0; i < _h_R_Durham[1]->size(); ++i) {
          IDataPoint* dp = _h_R_Durham[1]->point(i);
          double ycut = dp->coordinate(0)->value();
          if (y_23 < ycut && y_12 > ycut) {
            dp->coordinate(1)->setValue(dp->coordinate(1)->value() + weight);
          }
        }
        for (int i = 0; i < _h_R_Durham[2]->size(); ++i) {
          IDataPoint* dp = _h_R_Durham[2]->point(i);
          double ycut = dp->coordinate(0)->value();
          if (y_34 < ycut && y_23 > ycut) {
            dp->coordinate(1)->setValue(dp->coordinate(1)->value() + weight);
          }
        }
        for (int i = 0; i < _h_R_Durham[3]->size(); ++i) {
          IDataPoint* dp = _h_R_Durham[3]->point(i);
          double ycut = dp->coordinate(0)->value();
          if (y_45 < ycut && y_34 > ycut) {
            dp->coordinate(1)->setValue(dp->coordinate(1)->value() + weight);
          }
        }
        for (int i = 0; i < _h_R_Durham[4]->size(); ++i) {
          IDataPoint* dp = _h_R_Durham[4]->point(i);
          double ycut = dp->coordinate(0)->value();
          if (y_56 < ycut && y_45 > ycut) {
            dp->coordinate(1)->setValue(dp->coordinate(1)->value() + weight);
          }
        }
        for (int i = 0; i < _h_R_Durham[5]->size(); ++i) {
          IDataPoint* dp = _h_R_Durham[5]->point(i);
          double ycut = dp->coordinate(0)->value();
          if (y_56 > ycut) {
            dp->coordinate(1)->setValue(dp->coordinate(1)->value() + weight);
          }
        }
      }
    }


    void finalize() {
      if (!_initialised) return;
      
      normalize(_h_thrust);
      normalize(_h_heavyjetmass);
      normalize(_h_totaljetbroadening);
      normalize(_h_widejetbroadening);
      normalize(_h_cparameter);
      normalize(_h_thrustmajor);
      normalize(_h_thrustminor);
      normalize(_h_jetmassdifference);
      normalize(_h_aplanarity);
      normalize(_h_oblateness);
      normalize(_h_sphericity);
      
      
      for (size_t n = 0; n < 5; ++n) {
        scale(_h_y_Durham[n], 1.0/sumOfWeights());
      }
      
      for (size_t n = 0; n < 6; ++n) {
        // Scale integrated jet rates to 1
        for (int i = 0; i < _h_R_Durham[n]->size(); ++i) {
          IDataPoint* dp = _h_R_Durham[n]->point(i);
          dp->coordinate(1)->setValue(dp->coordinate(1)->value()/sumOfWeights());
        }
      }
      
    }


  private:
    
    bool _initialised;

    AIDA::IHistogram1D *_h_thrust;
    AIDA::IHistogram1D *_h_heavyjetmass;
    AIDA::IHistogram1D *_h_totaljetbroadening;
    AIDA::IHistogram1D *_h_widejetbroadening;
    AIDA::IHistogram1D *_h_cparameter;
    AIDA::IHistogram1D *_h_thrustmajor;
    AIDA::IHistogram1D *_h_thrustminor;
    AIDA::IHistogram1D *_h_jetmassdifference;
    AIDA::IHistogram1D *_h_aplanarity;
    AIDA::IHistogram1D *_h_oblateness;
    AIDA::IHistogram1D *_h_sphericity;
    
    AIDA::IDataPointSet *_h_R_Durham[6];
    AIDA::IHistogram1D *_h_y_Durham[5];

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<ALEPH_2004_S5765862> plugin_ALEPH_2004_S5765862;


}
