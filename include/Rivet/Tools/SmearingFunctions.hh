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


  //////////////////////////


  inline double TAU_EFF_ATLAS_RUN1(const Particle& t) {
    if (t.abseta() > 2.5) return 0;
    if (t.pT() < 1*GeV) return 0;
    /// @todo If leptonic return 0
    return 0.45;
  }


  /// From ATL-PHYS-PUB-2015-045 medium working point
  inline double TAU_EFF_ATLAS_RUN2(const Particle& t) {
    if (t.abseta() > 2.5) return 0;
    if (t.pT() < 1*GeV) return 0;
    /// @todo If leptonic return 0
    ///< @todo Make nProng-specific
    // 1 prong efficiency / mistag = 0.7 / 0.02
    // >= 2 prong efficiency / mistag = 0.6 / 0.01
    const vector<Particle> children = t.children();
    /// @todo Require prongs in tracker acceptance, too?
    const size_t nprongs = count_if(children.begin(), children.end(),
                                    [](const Particle& p){ return p.isHadron() && p.charge3() != 0; });
    if (nprongs == 0) return 0;
    if (nprongs == 1) return (t.abspid() == PID::TAU) ? 0.7 : 0.02;
    return (t.abspid() == PID::TAU) ? 0.6 : 0.01;
  }


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


  ///////////////////////


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

  //@}


}

#endif
