// -*- C++ -*-
#ifndef RIVET_CDF_2004_S5839831_HH
#define RIVET_CDF_2004_S5839831_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"


namespace Rivet {

  class CDF_2004_S5839831 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor: cuts on charged final state are \f$ -1 < \eta < 1 \f$ 
    /// and \f$ p_T > 0.4 \f$ GeV.
    CDF_2004_S5839831() 
      : _fsproj(-1.0, 1.0, 0.4*GeV), _jetproj(_fsproj, FastJets::CDFJETCLU, 0.7)
    {
      setBeams(PROTON, ANTIPROTON);
      
      addProjection(_fsproj);
      addProjection(_jetproj);
      
      /// @todo Declare that this is to be run on minimum bias data and jet 
      /// data with several ET triggers:
      /// * 1800 GeV: ET > 20, 50, 70 & 100 GeV
      /// * 630 GeV: ET > 5, 15 GeV
      /// Lots of runs needed to fill this paper!
    }

    /// Factory method
    static Analysis* create() { return new CDF_2004_S5839831(); }
    //@}

  public:

    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "5839831";
    }
    /// Get a description of the analysis.
    string getDescription() const {
      return "Transverse cone and 'Swiss cheese' CDF Run II underlying event analysis.";
    }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "CDF";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
      return "2004";
    }
    /// Journal, and preprint references.
    vector<string> getReferences() const {
      vector<string> ret;
      ret.push_back("Phys. Rev. D70, 072002 (2004)");
      ret.push_back("hep-ex/0404004");
      return ret;
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

    /// @name Internal projections
    //@{
    /// Use only charged tracks
    ChargedFinalState _fsproj;
    /// Lose 8% of charged tracks randomly.
    //LossyFinalState _fsproj;
    /// The jet algorithm used by this analysis.
    FastJets _jetproj;
    //@}

  private:

    /// @name Histogram collections
    //@{
    /// Profile histograms, binned in the \f$ E_T \f$ of the leading jet, for
    /// the \f$ p_T \f$ sum in the toward, transverse and away regions at 
    /// \f$ \sqrt{s} = 1800 \text{GeV} \f$.
    /// Corresponds to figure 2/3, and HepData table 2.
    AIDA::IProfile1D *_pt90Max1800,  *_pt90Min1800,  *_pt90Diff1800;

    /// Profile histograms, binned in the \f$ E_T \f$ of the leading jet, for
    /// the \f$ p_T \f$ sum in the toward, transverse and away regions at
    /// at \f$ \sqrt{s} = 630 \text{GeV} \f$.
    /// Corresponds to figure 8, and HepData table 8.
    AIDA::IProfile1D *_pt90Max630,   *_pt90Min630,   *_pt90Diff630;

    /// Profile histograms, binned in the \f$ E_T \f$ of the leading jet, for
    /// the cone track multiplicity at \f$ \sqrt{s} = 1800 \text{GeV} \f$.
    /// Corresponds to figure 5, and HepData table 4.
    AIDA::IProfile1D *_num90Max1800,  *_num90Min1800;

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
    AIDA::IHistogram1D *_numTracksDbn1800, *_ptDbn1800;
    //@}


  private:

    /// Hide the assignment operator.
    CDF_2004_S5839831& operator=(const CDF_2004_S5839831&);

  };

}

#endif
