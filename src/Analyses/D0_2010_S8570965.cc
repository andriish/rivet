// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {

  typedef std::pair<double, double> doublepair;

  class D0_2010_S8570965 : public Analysis {
  public:

    D0_2010_S8570965()
      : Analysis("D0_2010_S8570965") 
    {
      setBeams(PROTON, ANTIPROTON);
      setNeedsCrossSection(true);
    }


  public:

    void init() {
      FinalState fs;
      addProjection(fs, "FS");
   
      IdentifiedFinalState ifs(-0.9, 0.9, 20.0*GeV);
      ifs.acceptId(PHOTON);
      addProjection(ifs, "IFS");

      _h_M = bookHistogram1D(1, 1, 1);
      _h_pT = bookHistogram1D(2, 1, 1);
      _h_dPhi = bookHistogram1D(3, 1, 1);
      _h_costheta = bookHistogram1D(4, 1, 1);
      
      std::vector<doublepair> M_ranges;
      M_ranges.push_back(std::make_pair(30.0, 50.0));
      M_ranges.push_back(std::make_pair(50.0, 80.0));
      M_ranges.push_back(std::make_pair(80.0, 350.0));
      int i=0;
      
      foreach (const doublepair& M_range, M_ranges) {
        _h_pT_M.addHistogram(M_range.first, M_range.second, bookHistogram1D(5+i, 1, 1));
        _h_dPhi_M.addHistogram(M_range.first, M_range.second, bookHistogram1D(6+i, 1, 1));
        _h_costheta_M.addHistogram(M_range.first, M_range.second, bookHistogram1D(7+i, 1, 1));
        ++i;
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      ParticleVector photons = applyProjection<IdentifiedFinalState>(event, "IFS").particlesByPt();
      if (photons.size() < 2 ||
          (photons[0].momentum().pT() < 21.0*GeV)) {
        vetoEvent;
      }
      
      // Isolate photons with ET_sum in cone
      ParticleVector isolated_photons;
      ParticleVector fs = applyProjection<FinalState>(event, "FS").particles();
      foreach (const Particle& photon, photons) {
        FourMomentum mom_in_cone;
        double eta_P = photon.momentum().eta();
        double phi_P = photon.momentum().phi();
        double Etsum=0.0;
        foreach (const Particle& p, fs) {
          if (deltaR(eta_P, phi_P, p.momentum().eta(), p.momentum().phi()) < 0.4) {
            mom_in_cone += p.momentum();
            if (PID::threeCharge(p.pdgId())!=0) Etsum += p.momentum().Et();
          }
        }
        if (mom_in_cone.E()/photon.momentum().E() < 1.1 && Etsum<1.5*GeV) {
          isolated_photons.push_back(photon);
        }
      }
   
      if (isolated_photons.size() != 2) {
        vetoEvent;
      }
      
      FourMomentum y1=isolated_photons[0].momentum();
      FourMomentum y2=isolated_photons[0].momentum();
      if (deltaR(y1, y2)<0.4) {
        vetoEvent;
      }
      
      FourMomentum yy=y1+y2;
      double Myy = yy.mass()/GeV;
      double pTyy = yy.pT()/GeV;
      if (Myy<pTyy) {
        vetoEvent;
      }
      
      double dPhiyy = mapAngle0ToPi(y1.phi()-y2.phi());
      double costhetayy = fabs(tanh(y1.eta()-y2.eta())/2.0);
      
      _h_M->fill(Myy, weight);
      _h_pT->fill(pTyy, weight);
      _h_dPhi->fill(dPhiyy, weight);
      _h_costheta->fill(costhetayy, weight);
      
      _h_pT_M.fill(Myy, pTyy, weight);
      _h_dPhi_M.fill(Myy, dPhiyy, weight);
      _h_costheta_M.fill(Myy, costhetayy, weight);
    }


    void finalize() {

      scale(_h_M, crossSection()/sumOfWeights());
      scale(_h_pT, crossSection()/sumOfWeights());
      scale(_h_dPhi, crossSection()/sumOfWeights());
      scale(_h_costheta, crossSection()/sumOfWeights());
      std::vector<double> dMyy;
      dMyy += 20.0, 30.0, 270.0;
      for (size_t i=0; i<3; ++i) {
        scale(_h_pT_M.getHistograms()[i], crossSection()/sumOfWeights()/dMyy[i]);
        scale(_h_dPhi_M.getHistograms()[i], crossSection()/sumOfWeights()/dMyy[i]);
        scale(_h_costheta_M.getHistograms()[i], crossSection()/sumOfWeights()/dMyy[i]);
      }
      
    }


  private:

    AIDA::IHistogram1D *_h_M;
    AIDA::IHistogram1D *_h_pT;
    AIDA::IHistogram1D *_h_dPhi;
    AIDA::IHistogram1D *_h_costheta;
    BinnedHistogram<double> _h_pT_M;
    BinnedHistogram<double> _h_dPhi_M;
    BinnedHistogram<double> _h_costheta_M;

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<D0_2010_S8570965> plugin_D0_2010_S8570965;


}
