// -*- C++ -*-

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// Rivet analysis class for ATLAS 2010 minimum bias analysis
  class ATLAS_2010_S8918562 : public Analysis {
  public:

    /// Helper for collectively filling Nch, pT, eta, and pT vs. Nch histograms
    void fillPtEtaNch(const ChargedFinalState& cfs, int nchcut, double weight,
                      Histo1DPtr h_nch, Histo1DPtr h_pt,
                      Histo1DPtr h_eta, Profile1DPtr h_ptnch = Profile1DPtr()) {

      // Get number of particles and skip if event fails cut
      const int nch = cfs.size();
      if (nch < nchcut) return;

      // Fill nch
      h_nch->fill(nch, weight);
      // Loop over particles, fill pT, eta and ptnch
      foreach (const Particle& p, cfs.particles()) {
        const double pt = p.pT();
        h_pt->fill(pt/GeV, weight/pt);
        h_eta->fill(p.eta(), weight);
        if (h_ptnch != 0) h_ptnch->fill(nch, pt/GeV, weight);
      }
    }


    /// Default constructor
    ATLAS_2010_S8918562() : Analysis("ATLAS_2010_S8918562") {
      _sumW_pt100_nch2 = 0;
      _sumW_pt100_nch20 = 0;
      _sumW_pt500_nch1 = 0;
      _sumW_pt500_nch6 = 0;
      _sumW_pt2500_nch1 = 0;
    }


    /// Initialization, called once before running
    void init() {
      // Projections
      const ChargedFinalState cfs100(-2.5, 2.5, 100.0*MeV);
      declare(cfs100, "CFS100");
      const ChargedFinalState cfs500(-2.5, 2.5, 500.0*MeV);
      declare(cfs500, "CFS500");
      const ChargedFinalState cfs2500(-2.5, 2.5, 2500.0*MeV);
      declare(cfs2500, "CFS2500");

      // Book histograms
      if (fuzzyEquals(sqrtS()/GeV, 900)) {
        book(_hist_pt100_nch2_nch ,18, 1, 1);
        book(_hist_pt100_nch2_pt ,11, 1, 1);
        book(_hist_pt100_nch2_eta ,4, 1, 1);
        book(_hist_pt100_nch2_ptnch ,24, 1, 1);

        book(_hist_pt100_nch20_nch ,34, 1, 1);
        book(_hist_pt100_nch20_pt ,30, 1, 1);
        book(_hist_pt100_nch20_eta ,26, 1, 1);

        book(_hist_pt500_nch1_nch ,15, 1, 1);
        book(_hist_pt500_nch1_pt ,8, 1, 1);
        book(_hist_pt500_nch1_eta ,1, 1, 1);
        book(_hist_pt500_nch1_ptnch ,22, 1, 1);

        book(_hist_pt500_nch6_nch ,20, 1, 1);
        book(_hist_pt500_nch6_pt ,13, 1, 1);
        book(_hist_pt500_nch6_eta ,6, 1, 1);

        book(_hist_pt2500_nch1_nch ,36, 1, 1);
        book(_hist_pt2500_nch1_pt ,32, 1, 1);
        book(_hist_pt2500_nch1_eta ,28, 1, 1);
        book(_hist_pt2500_nch1_ptnch ,38, 1, 1);

      } else if (fuzzyEquals(sqrtS()/GeV, 2360)) {

        book(_hist_pt500_nch1_nch ,16, 1, 1);
        book(_hist_pt500_nch1_pt ,9, 1, 1);
        book(_hist_pt500_nch1_eta ,2, 1, 1);
        // This one histogram might be called while unbooked, so ensure its pointer is null!
        _hist_pt500_nch1_ptnch = Profile1DPtr();

      } else if (fuzzyEquals(sqrtS()/GeV, 7000)) {

        book(_hist_pt100_nch2_nch ,19, 1, 1);
        book(_hist_pt100_nch2_pt ,12, 1, 1);
        book(_hist_pt100_nch2_eta ,5, 1, 1);
        book(_hist_pt100_nch2_ptnch ,25, 1, 1);

        book(_hist_pt100_nch20_nch ,35, 1, 1);
        book(_hist_pt100_nch20_pt ,31, 1, 1);
        book(_hist_pt100_nch20_eta ,27, 1, 1);

        book(_hist_pt500_nch1_nch ,17, 1, 1);
        book(_hist_pt500_nch1_pt ,10, 1, 1);
        book(_hist_pt500_nch1_eta ,3, 1, 1);
        book(_hist_pt500_nch1_ptnch ,23, 1, 1);

        book(_hist_pt500_nch6_nch ,21, 1, 1);
        book(_hist_pt500_nch6_pt ,14, 1, 1);
        book(_hist_pt500_nch6_eta ,7, 1, 1);

        book(_hist_pt2500_nch1_nch ,37, 1, 1);
        book(_hist_pt2500_nch1_pt ,33, 1, 1);
        book(_hist_pt2500_nch1_eta ,29, 1, 1);
        book(_hist_pt2500_nch1_ptnch ,39, 1, 1);

      } else {
        throw LogicError("The ATLAS_2010_S8918562 analysis is only valid for sqrt(s) = 900, 2360 and 7000 GeV!");
      }

    }



    void analyze(const Event& event) {
      const double weight = 1.0;

      // 100 GeV final states
      if (!fuzzyEquals(sqrtS()/GeV, 2360)) {
        const ChargedFinalState& cfs100 = apply<ChargedFinalState>(event, "CFS100");
        // nch>=2
        if (cfs100.size() >= 2) _sumW_pt100_nch2 += weight;
        fillPtEtaNch(cfs100, 2, weight, _hist_pt100_nch2_nch, _hist_pt100_nch2_pt, _hist_pt100_nch2_eta, _hist_pt100_nch2_ptnch);
        // nch>=20
        if (cfs100.size() >= 20) _sumW_pt100_nch20 += weight;
        fillPtEtaNch(cfs100, 20, weight, _hist_pt100_nch20_nch, _hist_pt100_nch20_pt, _hist_pt100_nch20_eta);
      }

      // 500 GeV final states
      const ChargedFinalState& cfs500 = apply<ChargedFinalState>(event, "CFS500");
      // nch>=1
      if (cfs500.size() >= 1) _sumW_pt500_nch1 += weight;
      fillPtEtaNch(cfs500, 1, weight, _hist_pt500_nch1_nch, _hist_pt500_nch1_pt, _hist_pt500_nch1_eta, _hist_pt500_nch1_ptnch);
      // nch>=6
      if (!fuzzyEquals(sqrtS()/GeV, 2360)) {
        if (cfs500.size() >= 6) _sumW_pt500_nch6 += weight;
        fillPtEtaNch(cfs500, 6, weight, _hist_pt500_nch6_nch, _hist_pt500_nch6_pt, _hist_pt500_nch6_eta);
      }

      // 2500 GeV final states
      if (!fuzzyEquals(sqrtS()/GeV, 2360)) {
        const ChargedFinalState& cfs2500 = apply<ChargedFinalState>(event, "CFS2500");
        // nch>=1
        if (cfs2500.size() >= 1) _sumW_pt2500_nch1 += weight;
        fillPtEtaNch(cfs2500, 1, weight, _hist_pt2500_nch1_nch, _hist_pt2500_nch1_pt, _hist_pt2500_nch1_eta, _hist_pt2500_nch1_ptnch);
      }

    }



    void finalize() {

      if (!fuzzyEquals(sqrtS()/GeV, 2360)) {
        if (_sumW_pt100_nch2 > 0) {
          scale(_hist_pt100_nch2_nch, 1.0/_sumW_pt100_nch2);
          scale(_hist_pt100_nch2_pt, 1.0/_sumW_pt100_nch2/TWOPI/5);
          scale(_hist_pt100_nch2_eta, 1.0/_sumW_pt100_nch2);
        }

        if (_sumW_pt100_nch20 > 0) {
          scale(_hist_pt100_nch20_nch, 1.0/_sumW_pt100_nch20);
          scale(_hist_pt100_nch20_pt, 1.0/_sumW_pt100_nch20/TWOPI/5);
          scale(_hist_pt100_nch20_eta, 1.0/_sumW_pt100_nch20);
        }

        if (_sumW_pt500_nch6 > 0) {
          scale(_hist_pt500_nch6_nch, 1.0/_sumW_pt500_nch6);
          scale(_hist_pt500_nch6_pt, 1.0/_sumW_pt500_nch6/TWOPI/5);
          scale(_hist_pt500_nch6_eta, 1.0/_sumW_pt500_nch6);
        }

        if (_sumW_pt2500_nch1 > 0) {
          scale(_hist_pt2500_nch1_nch, 1.0/_sumW_pt2500_nch1);
          scale(_hist_pt2500_nch1_pt, 1.0/_sumW_pt2500_nch1/TWOPI/5);
          scale(_hist_pt2500_nch1_eta, 1.0/_sumW_pt2500_nch1);
        }
      }

      if (_sumW_pt500_nch1 > 0) {
        scale(_hist_pt500_nch1_nch, 1.0/_sumW_pt500_nch1);
        scale(_hist_pt500_nch1_pt, 1.0/_sumW_pt500_nch1/TWOPI/5);
        scale(_hist_pt500_nch1_eta, 1.0/_sumW_pt500_nch1);
      }
    }


  private:

    double _sumW_pt100_nch2, _sumW_pt100_nch20, _sumW_pt500_nch1,
      _sumW_pt500_nch6, _sumW_pt2500_nch1;

    Histo1DPtr _hist_pt100_nch2_nch,
      _hist_pt100_nch20_nch, _hist_pt500_nch1_nch,
      _hist_pt500_nch6_nch, _hist_pt2500_nch1_nch;

    Histo1DPtr _hist_pt100_nch2_pt,
      _hist_pt100_nch20_pt, _hist_pt500_nch1_pt,
      _hist_pt500_nch6_pt, _hist_pt2500_nch1_pt;

    Histo1DPtr _hist_pt100_nch2_eta,
      _hist_pt100_nch20_eta, _hist_pt500_nch1_eta,
      _hist_pt500_nch6_eta, _hist_pt2500_nch1_eta;

    Profile1DPtr _hist_pt100_nch2_ptnch,
      _hist_pt500_nch1_ptnch, _hist_pt2500_nch1_ptnch;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2010_S8918562);

}
