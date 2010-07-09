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


  /// @brief ALEPH jet rates and event shapes at LEP 1 and 2
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
      if (offset<6) { // there is no y34, y45 and y56 for 200 gev
        _h_y_Durham[2] = bookHistogram1D(offset+165, 1, 1); // y34 d165 ... d172, but not 171
        _h_y_Durham[3] = bookHistogram1D(offset+173, 1, 1); // y45 d173 ... d179
        _h_y_Durham[4] = bookHistogram1D(offset+180, 1, 1); // y56 d180 ... d186
      }
      else if (offset==6) {
        _h_y_Durham[2] = NULL;
        _h_y_Durham[3] = NULL;
        _h_y_Durham[4] = NULL;
      }
      else if (offset==7) {
        _h_y_Durham[2] = bookHistogram1D(172, 1, 1);
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
        for (size_t i=0; i<5; ++i) {
          if (_h_y_Durham[i]) {
            _h_y_Durham[i]->fill(-log(durjet.clusterSeq()->exclusive_ymerge_max(i+1)), weight);
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

      for (int N=1; N<7; ++N) {
        // calculate the N jet fraction from the jet resolution histograms

        for (size_t i = 0; i < _h_R_Durham[N-1]->size(); ++i) {
          IDataPoint* dp = _h_R_Durham[N-1]->point(i);
          // get ycut at which the njet-fraction is to be calculated
          /// @todo HepData has binwidths here, which doesn't make sense at all
          /// I assume the low edge to be the ycut
          double ycut = dp->coordinate(0)->value()-dp->coordinate(0)->errorMinus();

          // sum all >=N jet events
          double sigmaNinclusive = 0.0;
          if (N>1) {
            if (_h_y_Durham[N-2]) {
              AIDA::IHistogram1D* y_Nminus1_N = _h_y_Durham[N-2];
              // watch out, y_NM is negatively binned
              int cutbin=y_Nminus1_N->coordToIndex(-ycut);
              if (cutbin==AIDA::IAxis::UNDERFLOW_BIN) cutbin=0;
              if (cutbin==AIDA::IAxis::OVERFLOW_BIN) cutbin=y_Nminus1_N->axis().bins()-1;
              for (int ibin=0; ibin<cutbin; ++ibin) {
                sigmaNinclusive += y_Nminus1_N->binHeight(ibin);
              }
            }
          }
          else sigmaNinclusive = sumOfWeights();

          // sum all >=N+1 jet events
          double sigmaNplus1inclusive = 0.0;
          if (N<6) {
            if (_h_y_Durham[N-1]) {
              AIDA::IHistogram1D* y_N_Nplus1 = _h_y_Durham[N-1];
              // watch out, y_NM is negatively binned
              int cutbin=y_N_Nplus1->coordToIndex(-ycut);
              if (cutbin==AIDA::IAxis::UNDERFLOW_BIN) cutbin=0;
              if (cutbin==AIDA::IAxis::OVERFLOW_BIN) cutbin=y_N_Nplus1->axis().bins();
              for (int ibin=0; ibin<cutbin; ++ibin) {
                sigmaNplus1inclusive += y_N_Nplus1->binHeight(ibin);
              }
            }
          }

          // njetfraction = (sigma(>=N jet) - sigma(>=N+1 jet)) / sigma(tot)
          double njetfraction = (sigmaNinclusive-sigmaNplus1inclusive)/sumOfWeights();
          dp->coordinate(1)->setValue(njetfraction);
        }
      }


      for (size_t n = 0; n < 5; ++n) {
        if (_h_y_Durham[n]) {
          scale(_h_y_Durham[n], 1.0/sumOfWeights());
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
