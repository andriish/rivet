// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief ZEUS dijet photoproduction study used in the ZEUS jets PDF fit
  ///
  /// This class is a reproduction of the HZTool routine for the ZEUS
  /// dijet photoproduction paper which was used in the ZEUS jets PDF fit.
  ///
  /// @author Jon Butterworth
  class ZEUS_2001_S4815815 : public Analysis {
  public:

    /// Constructor
    ZEUS_2001_S4815815()
      : Analysis("ZEUS_2001_S4815815")
    {    }


    /// @name Analysis methods
    //@{

    // Book projections and histograms
    void init() {
      /// @todo Force conventional event rotation with proton along +z?
      FinalState fs;
      addProjection(FastJets(fs, FastJets::KT, 0.7), "Jets");

      IdentifiedFinalState positrons(fs, PID::POSITRON);
      addProjection(positrons, "Positrons");

      // Table 1
      _h_costh[0] = bookHisto1D(1, 1, 1);
      _h_costh[1] = bookHisto1D(1, 1, 2);
      // Table 2
      _h_etjet1[1][0] = bookHisto1D(2, 1, 1);
      _h_etjet1[1][1] = bookHisto1D(3, 1, 1);
      _h_etjet1[1][2] = bookHisto1D(4, 1, 1);
      _h_etjet1[1][3] = bookHisto1D(5, 1, 1);
      _h_etjet1[1][4] = bookHisto1D(6, 1, 1);
      _h_etjet1[1][5] = bookHisto1D(7, 1, 1);
      // Table 3
      _h_etjet1[0][0] = bookHisto1D(8, 1, 1);
      _h_etjet1[0][1] = bookHisto1D(9, 1, 1);
      _h_etjet1[0][2] = bookHisto1D(10, 1, 1);
      _h_etjet1[0][3] = bookHisto1D(11, 1, 1);
      _h_etjet1[0][4] = bookHisto1D(12, 1, 1);
      _h_etjet1[0][5] = bookHisto1D(13, 1, 1);
      // Table 4
      _h_etajet2[1][0] = bookHisto1D(14, 1, 1);
      _h_etajet2[1][1] = bookHisto1D(15, 1, 1);
      _h_etajet2[1][2] = bookHisto1D(16, 1, 1);
      // Table 5
      _h_etajet2[0][0] = bookHisto1D(17, 1, 1);
      _h_etajet2[0][1] = bookHisto1D(18, 1, 1);
      _h_etajet2[0][2] = bookHisto1D(19, 1, 1);
      // Table 6
      _h_xobsy[0] = bookHisto1D(20, 1, 1);
      _h_xobsy[1] = bookHisto1D(21, 1, 1);
      _h_xobsy[2] = bookHisto1D(22, 1, 1);
      _h_xobsy[3] = bookHisto1D(23, 1, 1);
    }


    // Do the analysis
    void analyze(const Event& event) {

      // Jet selection
      const Jets jets = applyProjection<FastJets>(event, "Jets").jets(Cuts::pT > 11*GeV && Cuts::etaIn(-1, 2.4), cmpMomByEt);
      if (jets.size() < 2) vetoEvent;
      if (jets[0].pT() < 14*GeV) vetoEvent;
      MSG_DEBUG("Jet multiplicity = " << jets.size());

      /// @todo Cut on pT/sqrt(Et)?

      // Leading jets and etabar & cos(theta*) computation
      const Jet& j1 = jets[0];
      const Jet& j2 = jets[1];
      const double etabar = (j1.eta() - j2.eta())/2.;
      const double costhetastar = tanh(etabar);

      // Get the scattered positron candidate
      const Particles positrons = applyProjection<FinalState>(event, "Positrons").particlesByPt();
      if (positrons.empty()) vetoEvent;
      const Particle& positron = positrons[0];
      const double Eeprime = positron.E();

      // Get the beams, for kinematics extraction
      const ParticlePair bs = beams();
      const double Ee = (bs.first.pid() == PID::POSITRON ? bs.first : bs.second).E();

      // Cuts on inelasticity y, and computation of x_y^obs
      const double inelasticity = 1 - (Eeprime/Ee/2.0)*(1 - positron.theta());
      if (!inRange(inelasticity, 0.2, 0.85)) vetoEvent; ///< @note The cuts on y_e and Y_JB look contradictory! These are just y_JB
      const double xyobs = (j1.Et() * exp(-j1.eta()) + j2.Et() * exp(-j2.eta())) / (2*inelasticity*Ee);
      const size_t i_xyobs = (xyobs < 0.75) ? 0 : 1;

      // Fill histograms
      const double weight = event.weight();
      // T1
      if ((j1.mom()+j2.mom()).mass() > 42*GeV && inRange(etabar, 0.1, 0.3))
        _h_costh[i_xyobs]->fill(abs(costhetastar), weight);
      // T2, T3
      if (inRange(j1.eta(), -1, 0) && inRange(j2.eta(), -1, 0))
        _h_etjet1[i_xyobs][0]->fill(j1.Et()/GeV, weight);
      else if (inRange(j1.eta(), 0, 1) && inRange(j2.eta(), -1, 0))
        _h_etjet1[i_xyobs][1]->fill(j1.Et()/GeV, weight);
      else if (inRange(j1.eta(), 0, 1) && inRange(j2.eta(), 0, 1))
        _h_etjet1[i_xyobs][2]->fill(j1.Et()/GeV, weight);
      else if (inRange(j1.eta(), 1, 2.4) && inRange(j2.eta(), -1, 0))
        _h_etjet1[i_xyobs][3]->fill(j1.Et()/GeV, weight);
      else if (inRange(j1.eta(), 1, 2.4) && inRange(j2.eta(), 0, 1))
        _h_etjet1[i_xyobs][4]->fill(j1.Et()/GeV, weight);
      else if (inRange(j1.eta(), 1, 2.4) && inRange(j2.eta(), 1, 2.4))
        _h_etjet1[i_xyobs][5]->fill(j1.Et()/GeV, weight);
      // T4, T5
      if (inRange(j1.eta(), -1, 0))
        _h_etajet2[i_xyobs][0]->fill(j2.eta(), weight);
      else if (inRange(j1.eta(), 0, 1))
        _h_etajet2[i_xyobs][1]->fill(j2.eta(), weight);
      else if (inRange(j1.eta(), 1, 2.4))
        _h_etajet2[i_xyobs][2]->fill(j2.eta(), weight);
      // T6
      if (inRange(j1.Et()/GeV, 14, 17))
        _h_xobsy[0]->fill(xyobs, weight);
      else if (inRange(j1.Et()/GeV, 17, 25))
        _h_xobsy[1]->fill(xyobs, weight);
      else if (inRange(j1.Et()/GeV, 25, 35))
        _h_xobsy[2]->fill(xyobs, weight);
      else if (inRange(j1.Et()/GeV, 35, 90))
        _h_xobsy[3]->fill(xyobs, weight);
    }


    // Finalize
    void finalize() {
      const double sf = crossSection()/picobarn/sumOfWeights();
      for (size_t ix = 0; ix < 2; ++ix) {
        scale(_h_costh[ix], sf);
        for (auto& h : _h_etjet1[ix]) scale(h, sf);
        for (auto& h : _h_etajet2[ix]) scale(h, sf);
      }
      for (auto& h : _h_xobsy) scale(h, sf);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_costh[2], _h_etjet1[2][6], _h_etajet2[2][3], _h_xobsy[4];
    //@}

  };


  DECLARE_RIVET_PLUGIN(ZEUS_2001_S4815815);

}
