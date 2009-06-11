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

  public:

    /// @name Publication metadata
    //@{
    /// A short description of the analysis.
    string spiresId() const {
      return "5839831";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Transverse cone and 'Swiss cheese' underlying event studies";
    }
    string description() const {
      ostringstream os;
      os << "This analysis studies the underlying event via transverse cones of "
         << " $R = 0.7$ at 90 degrees in \\phi relative to the leading (highest "
         << "E) jet, at sqrt(s) = 630 and 1800 GeV. This is similar to the 2001 "
         << "CDF UE analysis, except that cones, rather than the whole central "
         << "\\eta range are used. The transverse cones are categorised as TransMIN "
         << "and TransMAX on an event-by-event basis, to give greater sensitivity "
         << "to the UE component."
         << "\n\n"
         << "'Swiss Cheese' distributions, where cones around the leading $n$ "
         << "jets are excluded from the distributions, are also included for "
         << "$n = 2, 3$."
         << "\n\n"
         << "This analysis is useful for constraining the energy evolution of "
         << "the underlying event, since it performs the same analyses at two "
         << "distinct CoM energies."
         << "\n\n"
         << "WARNING: this analysis is not currently considered valid for MC "
         << "tuning and validation studies due to ambiguities in the paper and "
         << "non-reproducability of the MC plots shown in the paper. The fit to "
         << "data is sufficiently poor that this analysis skews the overall "
         << "goodness of fit in tuning studies, and has to be excluded. If you "
         << "can help to improve this analysis and make it usable for validation "
         << "studies, please get in touch!";
      return os.str();
    }

    /// Type of events required by this analysis
    string runInfo() const {
      ostringstream os;
      os << "* Two different beam energies: sqrt(s) = 630 & 1800 GeV\n"
         << "* Event type: generic QCD events\n"
         << "* Several pTmin cutoffs are probably required to fill the profile"
         << "  histograms, e.g.\n\n"
         << "  * { 0 (min bias), 30, 90, 150 GeV } at 1800 GeV; and\n"
         << "  * { 0 (min bias), 20, 90, 150 GeV } at 630 GeV";
      return os.str();
    }

    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "CDF";
    }
    /// Collider on which the experiment ran.
    string collider() const {
      return "Tevatron Run 2";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "2004";
    }
    /// Journal, and preprint references.
    vector<string> references() const {
      vector<string> ret;
      ret += "Phys. Rev. D70, 072002 (2004)";
      ret += "arXiv:hep-ex/0404004";
      return ret;
    }
    /// Routine authors
    vector<string> authors() const {
      vector<string> ret;
      ret += "Andy Buckley <andy.buckley@durham.ac.uk>";
      return ret;
    }
    /// Validation status
    string status() const {
      return "UNVALIDATED";
    }
    //@}

  public:

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
