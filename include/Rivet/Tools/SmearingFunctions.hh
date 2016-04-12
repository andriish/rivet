// -*- C++ -*-
#ifndef RIVET_SmearingFunctions_HH
#define RIVET_SmearingFunctions_HH

#include "Rivet/Particle.hh"
#include "Rivet/Jet.hh"
#include <random>

namespace Rivet {


  /// Return a uniformly sampled random number between 0 and 1
  /// @todo Move to (math?)utils
  /// @todo Need to isolate random generators to a single thread
  inline double rand01() {
    //return rand() / (double)RAND_MAX;
    static random_device rd;
    static mt19937 gen(rd());
    return generate_canonical<double, 10>(gen);
  }



  /// @name Standard particle efficiency and smearing functions
  //@{

  inline double PARTICLE_FN0(const Particle& p) { return 0; }
  inline double PARTICLE_FN1(const Particle& p) { return 1; }

  inline double P4_FN0(const FourMomentum& p) { return 0; }
  inline double P4_FN1(const FourMomentum& p) { return 1; }

  inline Particle PARTICLE_SMEAR_IDENTITY(const Particle& p) { return p; }


  ////////////////////////


  inline double ELECTRON_TRKEFF_ATLAS_RUN1(const Particle& e) {
    if (e.abseta() > 2.5) return 0;
    if (e.pT() < 0.1*GeV) return 0;
    if (e.abseta() < 1.5) {
      if (e.pT() < 1*GeV) return 0.73;
      if (e.pT() < 100*GeV) return 0.95;
      return 0.99;
    } else {
      if (e.pT() < 1*GeV) return 0.50;
      if (e.pT() < 100*GeV) return 0.83;
      else return 0.90;
    }
  }

  inline double ELECTRON_EFF_ATLAS_RUN1(const Particle& e) {
    if (e.abseta() > 2.5) return 0;
    if (e.pT() < 10*GeV) return 0;
    return (e.abseta() < 1.5) ? 0.95 : 0.85;
  }

  inline Particle ELECTRON_SMEAR_ATLAS_RUN1(const Particle& e) {
    /// @todo Need to isolate random generators to a single thread
    static random_device rd;
    static mt19937 gen(rd());

    static const vector<double> edges_eta = {0., 2.5, 3., 5.};
    static const vector<double> edges_pt = {0., 0.1, 25.};
    static const vector<double> e2s = {0.000, 0.015, 0.005,
                                       0.005, 0.005, 0.005,
                                       0.107, 0.107, 0.107};
    static const vector<double> es = {0.00, 0.00, 0.05,
                                      0.05, 0.05, 0.05,
                                      2.08, 2.08, 2.08};
    static const vector<double> cs = {0.00, 0.00, 0.25,
                                      0.25, 0.25, 0.25,
                                      0.00, 0.00, 0.00};

    const int i_eta = binIndex(e.abseta(), edges_eta, true);
    const int i_pt = binIndex(e.pT()/GeV, edges_pt, true);
    const int i = i_eta*edges_pt.size() + i_pt;

    // Calculate absolute resolution in GeV
    const double c1 = sqr(e2s[i]), c2 = sqr(es[i]), c3 = sqr(cs[i]);
    const double resolution = sqrt(c1*e.E2() + c2*e.E() + c3) * GeV;

    /// @todo Extract to a smear_energy helper function
    /// @todo Also make smear_direction and smear_pt functions, and jet versions that also update/smear constituents
    normal_distribution<> d(e.E(), resolution);
    const double smeared_E = max(d(gen), e.mass()); //< can't let the energy go below the mass!
    return Particle(e.pid(), FourMomentum::mkEtaPhiME(e.eta(), e.phi(), e.mass(), smeared_E));
  }


  inline double ELECTRON_TRKEFF_CMS_RUN1(const Particle& e) {
    if (e.abseta() > 2.5) return 0;
    if (e.pT() < 0.1*GeV) return 0;
    if (e.abseta() < 1.5) {
      return (e.pT() < 1*GeV) ? 0.70 : 0.95;
    } else {
      return (e.pT() < 1*GeV) ? 0.60 : 0.85;
    }
  }

  inline double ELECTRON_EFF_CMS_RUN1(const Particle& e) {
    if (e.abseta() > 2.5) return 0;
    if (e.pT() < 10*GeV) return 0;
    return (e.abseta() < 1.5) ? 0.95 : 0.85;
  }


  /// @brief CMS electron energy smearing, preserving direction
  ///
  /// Calculate resolution
  /// for pT > 0.1 GeV, E resolution = |eta| < 0.5 -> sqrt(0.06^2 + pt^2 * 1.3e-3^2)
  ///                                  |eta| < 1.5 -> sqrt(0.10^2 + pt^2 * 1.7e-3^2)
  ///                                  |eta| < 2.5 -> sqrt(0.25^2 + pt^2 * 3.1e-3^2)
  inline Particle ELECTRON_SMEAR_CMS_RUN1(const Particle& e) {
    /// @todo Need to isolate random generators to a single thread
    static random_device rd;
    static mt19937 gen(rd());

    // Calculate absolute resolution in GeV from functional form
    double resolution = 0;
    const double abseta = e.abseta();
    if (e.pT() > 0.1*GeV && abseta < 2.5) { //< should be a given from efficiencies
      if (abseta < 0.5) {
        resolution = add_quad(0.06, 1.3e-3 * e.pT()/GeV) * GeV;
      } else if (abseta < 1.5) {
        resolution = add_quad(0.10, 1.7e-3 * e.pT()/GeV) * GeV;
      } else { // still |eta| < 2.5
        resolution = add_quad(0.25, 3.1e-3 * e.pT()/GeV) * GeV;
      }
    }

    /// @todo Extract to a smear_energy helper function
    normal_distribution<> d(e.E(), resolution);
    const double smeared_E = max(d(gen), e.mass()); //< can't let the energy go below the mass!
    return Particle(e.pid(), FourMomentum::mkEtaPhiME(e.eta(), e.phi(), e.mass(), smeared_E));
  }


  ///////////////////


  inline double MUON_TRKEFF_ATLAS_RUN1(const Particle& m) {
    if (m.abseta() > 2.5) return 0;
    if (m.pT() < 0.1*GeV) return 0;
    if (m.abseta() < 1.5) {
      return (m.pT() < 1*GeV) ? 0.75 : 0.99;
    } else {
      return (m.pT() < 1*GeV) ? 0.70 : 0.98;
    }
  }

  inline double MUON_EFF_ATLAS_RUN1(const Particle& m) {
    if (m.abseta() > 2.7) return 0;
    if (m.pT() < 10*GeV) return 0;
    return (m.abseta() < 1.5) ? 0.95 : 0.85;
  }

  inline Particle MUON_SMEAR_ATLAS_RUN1(const Particle& m) {
    /// @todo Need to isolate random generators to a single thread
    static random_device rd;
    static mt19937 gen(rd());

    static const vector<double> edges_eta = {0, 1.5, 2.5};
    static const vector<double> edges_pt = {0, 0.1, 1.0, 10., 200.};
    static const vector<double> res = {0., 0.03, 0.02, 0.03, 0.05,
                                       0., 0.04, 0.03, 0.04, 0.05};

    const int i_eta = binIndex(m.abseta(), edges_eta, true);
    const int i_pt = binIndex(m.pT()/GeV, edges_pt, true);
    const int i = i_eta*edges_pt.size() + i_pt;

    const double resolution = res[i];

    // Smear by a Gaussian centered on the current pT, with width given by the resolution
    /// @todo Extract to a smear_pt helper function
    normal_distribution<> d(m.pT(), resolution*m.pT());
    const double smeared_pt = max(d(gen), 0.);
    return Particle(m.pid(), FourMomentum::mkEtaPhiMPt(m.eta(), m.phi(), m.mass(), smeared_pt));
  }


  /// @note Eff values currently identical to those in ATLAS (AB, 2016-04-12)
  inline double MUON_TRKEFF_CMS_RUN1(const Particle& m) {
    if (m.abseta() > 2.5) return 0;
    if (m.pT() < 0.1*GeV) return 0;
    if (m.abseta() < 1.5) {
      return (m.pT() < 1*GeV) ? 0.75 : 0.99;
    } else {
      return (m.pT() < 1*GeV) ? 0.70 : 0.98;
    }
  }


  inline double MUON_EFF_CMS_RUN1(const Particle& m) {
    if (m.abseta() > 2.4) return 0;
    if (m.pT() < 10*GeV) return 0;
    return 0.95 * (m.abseta() < 1.5 ? 1 : exp(0.5 - 5e-4*m.pT()/GeV));
  }


  inline Particle MUON_SMEAR_CMS_RUN1(const Particle& m) {
    /// @todo Need to isolate random generators to a single thread
    static random_device rd;
    static mt19937 gen(rd());

    // Calculate fractional resolution
    // for pT > 0.1 GeV, mom resolution = |eta| < 0.5 -> sqrt(0.01^2 + pt^2 * 2.0e-4^2)
    //                                    |eta| < 1.5 -> sqrt(0.02^2 + pt^2 * 3.0e-4^2)
    //                                    |eta| < 2.5 -> sqrt(0.05^2 + pt^2 * 2.6e-4^2)
    double resolution = 0;
    const double abseta = m.abseta();
    if (m.pT() > 0.1*GeV && abseta < 2.5) {
      if (abseta < 0.5) {
        resolution = add_quad(0.01, 2.0e-4 * m.pT()/GeV);
      } else if (abseta < 1.5) {
        resolution = add_quad(0.02, 3.0e-4 * m.pT()/GeV);
      } else { // still |eta| < 2.5... but isn't CMS' mu acceptance < 2.4?
        resolution = add_quad(0.05, 2.6e-4 * m.pT()/GeV);
      }
    }

    // Smear by a Gaussian centered on the current pT, with width given by the resolution
    /// @todo Extract to a smear_pt helper function
    normal_distribution<> d(m.pT(), resolution*m.pT());
    const double smeared_pt = max(d(gen), 0.);
    return Particle(m.pid(), FourMomentum::mkEtaPhiMPt(m.eta(), m.phi(), m.mass(), smeared_pt));
  }


  //////////////////////////


  /// @brief ATLAS Run 1 8 TeV tau efficiencies (medium working point)
  ///
  /// Taken from http://arxiv.org/pdf/1412.7086.pdf
  ///   20-40 GeV 1-prong LMT eff|mis = 0.66|1/10, 0.56|1/20, 0.36|1/80
  ///   20-40 GeV 3-prong LMT eff|mis = 0.45|1/60, 0.38|1/100, 0.27|1/300
  ///   > 40 GeV 1-prong LMT eff|mis = 0.66|1/15, 0.56|1/25, 0.36|1/80
  ///   > 40 GeV 3-prong LMT eff|mis = 0.45|1/250, 0.38|1/400, 0.27|1/1300
  inline double TAU_EFF_ATLAS_RUN1(const Particle& t) {
    if (t.abseta() > 2.5) return 0; //< hmm... mostly
    double pThadvis = 0;
    Particles chargedhadrons;
    for (const Particle& p : t.children()) {
      if (p.isHadron()) {
        pThadvis += p.pT(); //< right definition? Paper is unclear
        if (p.charge3() != 0 && p.abseta() < 2.5 && p.pT() > 1*GeV) chargedhadrons += p;
      }
    }
    if (chargedhadrons.empty()) return 0; //< leptonic tau
    if (pThadvis < 20*GeV) return 0; //< below threshold
    if (pThadvis < 40*GeV) {
      if (chargedhadrons.size() == 1) return (t.abspid() == PID::TAU) ? 0.56 : 1/20.;
      if (chargedhadrons.size() == 3) return (t.abspid() == PID::TAU) ? 0.38 : 1/100.;
    } else {
      if (chargedhadrons.size() == 1) return (t.abspid() == PID::TAU) ? 0.56 : 1/25.;
      if (chargedhadrons.size() == 3) return (t.abspid() == PID::TAU) ? 0.38 : 1/400.;
    }
    return 0;
  }


  /// @brief ATLAS Run 2 13 TeV tau efficiencies (medium working point)
  ///
  /// From https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PUBNOTES/ATL-PHYS-PUB-2015-045/ATL-PHYS-PUB-2015-045.pdf
  ///   LMT 1 prong efficiency/mistag = 0.6|1/30, 0.55|1/50, 0.45|1/120
  ///   LMT 3 prong efficiency/mistag = 0.5|1/30, 0.4|1/110, 0.3|1/300
  inline double TAU_EFF_ATLAS_RUN2(const Particle& t) {
    if (t.abseta() > 2.5) return 0; //< hmm... mostly
    double pThadvis = 0;
    Particles chargedhadrons;
    for (const Particle& p : t.children()) {
      if (p.isHadron()) {
        pThadvis += p.pT(); //< right definition? Paper is unclear
        if (p.charge3() != 0 && p.abseta() < 2.5 && p.pT() > 1*GeV) chargedhadrons += p;
      }
    }
    if (chargedhadrons.empty()) return 0; //< leptonic tau
    if (pThadvis < 20*GeV) return 0; //< below threshold
    if (chargedhadrons.size() == 1) return (t.abspid() == PID::TAU) ? 0.55 : 1/50.;
    if (chargedhadrons.size() == 3) return (t.abspid() == PID::TAU) ? 0.40 : 1/110.;
    return 0;
  }


  /// ATLAS Run 1 tau smearing
  /// @todo Currently a copy of the crappy jet smearing that is probably wrong...
  inline Particle TAU_SMEAR_ATLAS_RUN1(const Particle& t) {
    /// @todo Need to isolate random generators to a single thread
    static random_device rd;
    static mt19937 gen(rd());

    // Const fractional resolution for now
    static const double resolution = 0.03;

    // Smear by a Gaussian centered on 1 with width given by the (fractional) resolution
    /// @todo Is this the best way to smear? Should we preserve the energy, or pT, or direction?
    normal_distribution<> d(1., resolution);
    const double fsmear = max(d(gen), 0.);
    return Particle(t.pid(), FourMomentum::mkXYZM(t.px()*fsmear, t.py()*fsmear, t.pz()*fsmear, t.mass()));
  }


  /// CMS Run 2 tau efficiency
  ///
  /// @todo Needs work; this is the dumb version from Delphes 3.3.2
  inline double TAU_EFF_CMS_RUN2(const Particle& t) {
    return (t.abspid() == PID::TAU) ? 0.6 : 0;
  }

  /// CMS Run 1 tau smearing
  /// @todo Currently a copy of the crappy ATLAS one
  inline Particle TAU_SMEAR_CMS_RUN1(const Particle& t) {
    return TAU_SMEAR_ATLAS_RUN1(t);
  }

  /// CMS Run 2 tau smearing
  /// @todo Currently a copy of the Run 1 version
  inline Particle TAU_SMEAR_CMS_RUN2(const Particle& t) {
    return TAU_SMEAR_CMS_RUN1(t);
  }

  //@}


  /// @name Standard jet efficiency and smearing functions
  //@{

  /// Return a constant 0 given a Jet as argument
  inline double JET_EFF_ZERO(const Jet& p) { return 0; }
  /// Return a constant 1 given a Jet as argument
  inline double JET_EFF_ONE(const Jet& p) { return 1; }

  /// Return 1 if the given Jet contains a b, otherwise 0
  inline double JET_BTAG_PERFECT(const Jet& j) { return j.bTagged() ? 1 : 0; }
  /// Return the ATLAS Run 1 jet flavour tagging efficiency for the given Jet
  inline double JET_BTAG_ATLAS_RUN1(const Jet& j) {
    if (j.bTagged()) return 0.80*tanh(0.003*j.pT()/GeV)*(30/(1+0.086*j.pT()/GeV));
    if (j.cTagged()) return 0.20*tanh(0.02*j.pT()/GeV)*(1/(1+0.0034*j.pT()/GeV));
    return 0.002 + 7.3e-6*j.pT()/GeV;
  }

  /// Return 1 if the given Jet contains a c, otherwise 0
  inline double JET_CTAG_PERFECT(const Jet& j) { return j.cTagged() ? 1 : 0; }

  /// Take a jet and return an unmodified copy
  /// @todo Modify constituent particle vectors for consistency
  /// @todo Set a null PseudoJet if the Jet is smeared?
  inline Jet JET_SMEAR_IDENTITY(const Jet& j) { return j; }

  /// ATLAS Run 1 jet smearing
  /// @todo This is a cluster-level flat 3% resolution, I think, and smearing is suboptimal: improve!
  inline Jet JET_SMEAR_ATLAS_RUN1(const Jet& j) {
    /// @todo Need to isolate random generators to a single thread
    static random_device rd;
    static mt19937 gen(rd());

    // Const fractional resolution for now
    static const double resolution = 0.03;

    // Smear by a Gaussian centered on 1 with width given by the (fractional) resolution
    /// @todo Is this the best way to smear? Should we preserve the energy, or pT, or direction?
    normal_distribution<> d(1., resolution);
    const double fsmear = max(d(gen), 0.);
    return Jet(FourMomentum::mkXYZM(j.px()*fsmear, j.py()*fsmear, j.pz()*fsmear, j.mass()));
  }

  /// CMS Run 1 jet smearing
  /// @todo Just a copy of the suboptimal ATLAS one: improve!!
  inline Jet JET_SMEAR_CMS_RUN1(const Jet& j) {
    return JET_SMEAR_ATLAS_RUN1(j);
  }

  //@}


}

#endif
