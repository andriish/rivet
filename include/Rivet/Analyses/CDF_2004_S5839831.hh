// -*- C++ -*-
#ifndef RIVET_CDF_2004_S5839831_HH
#define RIVET_CDF_2004_S5839831_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Jet.hh"

namespace Rivet {

  /**
   * @brief "Acosta" CDF underlying event analysis
   * @author Andy Buckley
   */
  class CDF_2004_S5839831 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor: cuts on charged final state are \f$ -1 < \eta < 1 \f$ 
    /// and \f$ p_T > 0.4 \f$ GeV.
    CDF_2004_S5839831();

    /// Factory method
    static Analysis* create() { 
      return new CDF_2004_S5839831(); 
    }
    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}


  private:

    struct ConesInfo {
      ConesInfo() : numMax(0), numMin(0), ptMax(0), ptMin(0), ptDiff(0) {}
      unsigned int numMax, numMin;
      double ptMax, ptMin, ptDiff;
    };

    const ConesInfo calcTransCones(const double etaLead, const double phiLead, const ParticleVector& tracks);
    const ConesInfo calcTransCones(const FourMomentum& leadvec, const ParticleVector& tracks);

  private:

    /// @name Histogram collections
    //@{
    /// Profile histograms, binned in the \f$ E_T \f$ of the leading jet, for
    /// the average \f$ p_T \f$ in the toward, transverse and away regions at 
    /// \f$ \sqrt{s} = 1800 \text{GeV} \f$.
    /// Corresponds to Table 1, and HepData table 1.
    AIDA::IProfile1D *_pt90MaxAvg1800, *_pt90MinAvg1800;

    /// Profile histograms, binned in the \f$ E_T \f$ of the leading jet, for
    /// the \f$ p_T \f$ sum in the toward, transverse and away regions at 
    /// \f$ \sqrt{s} = 1800 \text{GeV} \f$.
    /// Corresponds to figure 2/3, and HepData table 2.
    AIDA::IProfile1D *_pt90Max1800, *_pt90Min1800, *_pt90Diff1800;

    /// Profile histograms, binned in the \f$ E_T \f$ of the leading jet, for
    /// the \f$ p_T \f$ sum in the toward, transverse and away regions at
    /// at \f$ \sqrt{s} = 630 \text{GeV} \f$.
    /// Corresponds to figure 8, and HepData table 8.
    AIDA::IProfile1D *_pt90Max630, *_pt90Min630, *_pt90Diff630;

    /// Profile histograms, binned in the \f$ E_T \f$ of the leading jet, for
    /// the cone track multiplicity at \f$ \sqrt{s} = 1800 \text{GeV} \f$.
    /// Corresponds to figure 5, and HepData table 4.
    AIDA::IProfile1D *_num90Max1800, *_num90Min1800;

    /// Profile histograms, binned in the \f$ E_T \f$ of the leading jet, for
    /// the \f$ p_T \f$ sum at \f$ \sqrt{s} = 1800 \text{GeV} \f$.
    /// Corresponds to figure 7, and HepData table 7.
    AIDA::IProfile1D *_pTSum1800_2Jet, *_pTSum1800_3Jet;

    /// Profile histograms, binned in the \f$ E_T \f$ of the leading jet, for
    /// the \f$ p_T \f$ sum at \f$ \sqrt{s} = 630 \text{GeV} \f$.
    /// Corresponds to figure 9, and HepData table 9.
    AIDA::IProfile1D *_pTSum630_2Jet, *_pTSum630_3Jet;

    /// Histogram of \f$ p_{T\text{sum}} \f$ distribution for 5 different 
    /// \f$ E_{T1} \f$ bins.
    /// Corresponds to figure 4, and HepData table 3.
    AIDA::IHistogram1D *_pt90Dbn1800Et40, *_pt90Dbn1800Et80, *_pt90Dbn1800Et120, 
      *_pt90Dbn1800Et160, *_pt90Dbn1800Et200;

    /// Histograms of track multiplicity and \f$ p_T \f$ distributions for 
    /// minimum bias events.
    /// Figure 6, and HepData tables 5 & 6.
    /// Figure 10, and HepData tables 10 & 11.
    AIDA::IHistogram1D *_numTracksDbn1800MB, *_ptDbn1800MB;
    AIDA::IHistogram1D *_numTracksDbn630MB, *_ptDbn630MB;
    //@}


    //private:
    //UniformRealRNG _rngEtaMB, _rngPhiMB;

  };


}

#endif
