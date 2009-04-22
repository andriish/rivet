// -*- C++ -*-
#ifndef RIVET_CDF_2001_S4751469_HH
#define RIVET_CDF_2001_S4751469_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /* @brief "Field-Stuart" CDF Run I underlying event analysis
   * @author Andy Buckley
   * 
   * 
   * @par Run conditions
   * 
   * @arg \f$ \sqrt{s} = \f$ 1800 GeV
   * @arg Run with generic QCD events.
   * @arg Several \f$ p_\perp^\text{min} \f$ cutoffs are probably required to fill the profile histograms:
   *   @arg \f$ p_\perp^\text{min} = \f$ 0 (min bias), 10, 20 GeV
   * 
   */ 
  class CDF_2001_S4751469 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor: cuts on final state are \f$ -1 < \eta < 1 \f$ 
    /// and \f$ p_T > 0.5 \f$ GeV. Use a lossy charged FS projection, which
    /// randomly discards 8% of charged particles, as a kind of hacky detector 
    /// correction.
    CDF_2001_S4751469();

    /// Factory method
    static Analysis* create() {
      return new CDF_2001_S4751469();
    }
    //@}


  public:

    /// @name Publication metadata
    //@{
    /// Get the SPIRES ID
    string spiresId() const {
      return "4751469";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Field & Stuart Run I underlying event analysis.";
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "CDF";
    }
    /// Collider on which the experiment ran.
    string collider() const {
      return "Tevatron Run 1";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "2001";
    }
    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> rtn;
      rtn += "Andy Buckley <andy.buckley@durham.ac.uk>";
      return rtn;
    }
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "The original CDF underlying event analysis, based on decomposing each "
         << "event into a transverse structure with ``toward'', ``away'' and "
         << "``transverse'' regions defined relative to the azimuthal direction of "
         << "the leading jet in the event. Since the toward region is by definition "
         << "dominated by the hard process, as is the away region by momentum balance "
         << "in the matrix element, the transverse region is most sensitive to "
         << "multi-parton interactions. The transverse regions occupy "
         << "$|\\phi| \\in [60\\degree, 120\\degree]$ for $|\\eta| < 1$. The pT "
         << "ranges for the leading jet are divided experimentally into the `min-bias' "
         << "sample from 0--20 GeV, and the `JET20' sample from 18--49 GeV.";
      return os.str();
    }
    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "CDF Run I conditions: ppbar QCD interactions at 1800 GeV. "
         << "Leading jet bins from 0--49 GeV: usually can be filled with a single "
         << "generator run without kinematic sub-samples, with $\\sim 1M$ events.";
      return os.str();
    }
    string status() const {
      return "VALIDATED";
    }
    /// Journal, and preprint references.
    vector<string> references() const {
      vector<string> ret;
      ret += "Phys.Rev.D65:092002,2002", "FNAL-PUB 01/211-E";
      // No arXiv code
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

    /// Sum total number of charged particles in the trans region, in 3 \f$ p_\perp^\text{lead} \f$ bins.
    double _totalNumTrans2, _totalNumTrans5, _totalNumTrans30;

    /// Sum the total number of events in 3 \f$ p_\perp^\text{lead} \f$ bins.
    double _sumWeightsPtLead2,_sumWeightsPtLead5, _sumWeightsPtLead30;

  private:

    /// @name Histogram collections
    //@{
    /// Profile histograms, binned in the \f$ p_T \f$ of the leading jet, for
    /// the \f$ p_T \f$ sum in the toward, transverse and away regions.
    AIDA::IProfile1D *_ptsumTowardMB,  *_ptsumTransMB,  *_ptsumAwayMB;
    AIDA::IProfile1D *_ptsumTowardJ20, *_ptsumTransJ20, *_ptsumAwayJ20;

    /// Profile histograms, binned in the \f$ p_T \f$ of the leading jet, for
    /// the number of charged particles per jet in the toward, transverse and
    /// away regions.
    AIDA::IProfile1D *_numTowardMB,  *_numTransMB,  *_numAwayMB;
    AIDA::IProfile1D *_numTowardJ20, *_numTransJ20, *_numAwayJ20;

    /// Histogram of \f$ p_T \f$ distribution for 3 different \f$ p_{T1} \f$ IR cutoffs.
    AIDA::IHistogram1D *_ptTrans2, *_ptTrans5, *_ptTrans30;
    //@}

  };


}

#endif
