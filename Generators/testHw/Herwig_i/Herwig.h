// --------------------------------------------------
//
// File:  Generators/Herwig.h
// Description:
//    This code is used to get a HERWIG Monte Carlo event.
//    genInitialize() is used to read parameters
//    callGenerator() makes the event
//    genFinalize() writes log files etc
//    fillEvt(GeneratorEvent* evt) passes the event to HepMC
//
//    The output will be stored in the transient event store so it can be
//    passed to the simulation.
//
// AuthorList:
//         Georgios Stavropoulos January 2003 : Migrate to Herwig 6.5
//         Ian Hinchliffe Sept 2000
//         Modeled after the CDF code by Marge Shapiro, Elena Bucciol, and  j ly
// Borut Paul Kersevan
// added access to lesatlas common (not part of herwig but needed for external procees February 2003

#ifndef GENERATORMODULESHERWIG_H
#define GENERATORMODULESHERWIG_H

#include "HepMC/GenEvent.h"

//#include "GeneratorModules/GenModule.h"
//#include "AcerMC_i/AcerMC_acset.h"
//#include "Tauola_i/Atlas_HEPEVT.h"

#include <string>
#include <vector>

typedef std::vector<std::string> CommandVector;

class StatusCode {
public:
  StatusCode() { _value = 0; }
  StatusCode(int newval) { _value = newval; }
  StatusCode(const StatusCode& newsc) { _value = newsc.value(); }

  int value() const { return _value; }
  int operator()() const { return this->value(); }

  static const int SUCCESS = 0;
  static const int FAILURE = 1;

private:
  int _value;
};

//class Herwig : public GenModule {
class Herwig {
public:
  Herwig(const std::string& name);
  virtual ~Herwig();
  virtual StatusCode genInitialize();
  virtual StatusCode callGenerator();
  virtual StatusCode genFinalize();
  virtual StatusCode fillEvt(HepMC::GenEvent* evt);
  void Localpar2Her();
  void Her2Localpar();
  // accessor for external processes
  //*AGB* AcerMC_acset& acermc_acset();


private:
  // data members corresponding to Herwig variables
  int m_process;              // QCD process generated
  int m_maxer_Parm;           //
  int m_randomseed1_Parm;      //
  int m_randomseed2_Parm;      //
  //cjl  <float>m_lambdaQCD_Parm;     // to calculate alphas
  double m_lambdaQCD_Parm;     // to calculate alphas
  //  string m_decayTable_Parm;   // parm for the file name that has the decay table
  //  int  m_scaleChoiceParm;     // choice of QCD scale: 1=pt/2, 2=pt, 3=2pt

  // commands of Menu 01
  std::string m_Beamtype1_Parm;// initialize beam particle 1
  double m_Energy1_Parm;//initialize beam energy 1
  std::string m_Beamtype2_Parm;// initialize beam particle 2
  double m_Energy2_Parm;//initialize beam energy 2

  // commands of Menu 02 - v6.1 no read/write
  int m_sudord_Parm ; //Sudakov order

  // commands of Menu 03
  std::string m_autpdf;
  int m_modpdf;
  double m_pltcut_Parm;

  // commands of Menu 04
  //cjl  double m_rapiditiesm_Parm;// select rapidity limits min
  double m_ptm_Parm;// select transverse momentum limits min
  double m_invarmassm_Parm;// select invariant mass limits min
  double m_rapiditiesM_Parm;// select rapidity limits max
  double m_ptM_Parm;// select transverse momentum limits max
  double m_invarmassM_Parm;// select invariant mass limits max
  double m_thrust_Parm;// select thrust
  double m_ptpow_Parm;// select ET power
  double m_empow_Parm;// select M power
  double m_q2pow_Parm;// select Q**2 power
  double m_qlim_Parm;// select upper limits in scale
  double m_q2dilsm_Parm;// set Q**2 limits min
  double m_q2dilsM_Parm;// set Q**2 limits max
  int m_bgshat_Parm;// scaling switch

  // commands of Menu 05
  int m_azsoft_Parm;//azimuthal correlation soft gluon
  int m_prsof_Parm;
  int m_azspin_Parm;//azimuthal correlation spin
  double m_qspac_Parm;//lower limit for spacelike evolution
  int m_ispac_Parm;//controls exact meaning of qspac
  int m_nospac_Parm;//can turn off ISR
  double m_vqcut_Parm;//quark virtuality cutoff
  double m_vgcut_Parm;//gluon virtuality cutoff
  double m_vpcut_Parm;//photon virtuality cutoff
  //  double m_gluonmass_Parm;//gluon virtual mass cutoff
  double m_intrinsicpt_Parm;// intrinsic pt incoming hadrons

  // commands of Menu 06
  double m_clmax_Parm;//mass cluster mass parameter
  double m_clpow_Parm;// mass dependence of cluster hadronization
  double m_psplt1_Parm;// mass dist in cluster splitting, non-b
  double m_psplt2_Parm;// mass dist in cluster splitting, b
  int m_cldir1_Parm;// direction memory, cluster splitting, non-b
  int m_cldir2_Parm;// direction memory, cluster splitting, b
  double m_clsmr1_Parm;// alignment, cluster splitting, non-b
  double m_clsmr2_Parm;// alignment, cluster splitting, b
  double m_qdiqk_Parm; //scale at which gluon can be split into diquarks
  double m_pdiqk_Parm; //probability of diquark splitting per unit log scale
  double m_decwt_Parm; // a priori weights for decuplet baryons for cluster decays
  double m_pwt_Parm; // a priori weights for quarks and diquarks

  // commands of Menu 07
  double m_pmbn1_Parm;// soft int. param for multiplicity
  double m_pmbn2_Parm;// soft int. param for multiplicity
  double m_pmbn3_Parm;// soft int. param for multiplicity
  double m_pmbk1_Parm;// soft int. param for negative binomial
  double m_pmbk2_Parm;// soft int. param for negative binomial
  double m_pmbm1_Parm;// soft int. param for cluster mass
  double m_pmbm2_Parm;// soft int. param for cluster mass
  double m_pmbp1_Parm;// soft int. param for cluster pt
  double m_pmbp2_Parm;// soft int. param for cluster pt
  double m_pmbp3_Parm;// soft int. param for cluster pt
  double m_ensof_Parm;//multy. enhancement for underlyng soft event

  // commands of Menu 08
  int m_pizstable_Parm; //whether pizero is stable (non-decaying)

  // commands of Menu 09
  int m_eventweight_Parm; //control event weights
  double m_maxweight_Parm; //set max weight
  double m_effmin_Parm; //min. acceptable eff to generate evts
  bool m_noshower_Parm; //turn off parton shower - changed to bool ((IH)
  bool m_nohadroniz_Parm; //turn off hadronization -- changed to bool (IH)
  double m_tlout_Parm; //time (sec) needed to terminate gracefully

  // commands of Menu 10
  int m_colors_Parm; // number of colors
  int m_flavors_Parm; // number of flavors
  double m_weiangle_Parm; // sin**2 of Weimberg angle
  double m_cabangle_Parm; // sin**2 of Cabibbo angle
  double m_wwidth_Parm; // Width of W boson
  double m_zwidth_Parm; //Width of Z boson
  double m_hwidth_Parm; //Width of Higgs boson

  // commands of Menu 11
  double m_down_Parm; // Set down quark mass
  double m_up_Parm; // Set up quark mass
  double m_strange_Parm; // Set strange quark mass
  double m_charm_Parm; // Set charm quark mass
  double m_bottom_Parm; // Set bottom quark mass
  double m_top_Parm; // Set top quark mass
  double m_gluonmass_Parm; // Set gluon mass
  double m_W_Parm; // Set W+/- mass
  double m_Z0_Parm; // Set Z0  mass
  double m_higgs_Parm; // Set Higgs  mass
  int    m_modBos_Parm[9]; // SB: Set decay modes for Weak bosons

  // commands of Menu 12
  std::string m_read_Filesusy; //parm for the susy file name

  // commands of Menu 13
  int m_zprime_Parm; // Whether to use Z-prime particle
  double m_ZP_Parm; // Mass of Z-prime particle
  double m_zpwidth_Parm; // Width of Z-prime particle

  // commands of Menu 14 
  int m_iprint_Parm; // Herwig print level control
  int m_maxpr_Parm; // Dump of Herwig events  
  int m_prvtx_Parm; // Whether to have vtx print in dump  
  int m_prndef_Parm; // print ASCII to screen/log file  
  int m_prntex_Parm; // print LATEX to *.tex files  
  int m_prnweb_Parm; // print html to *.html  
  int m_nprfmt_Parm; // controls print sig. figs - short or long  
  //  
  // graviton stuff  
  double m_grvlam_Parm;  double m_emgrv_Parm;  double m_gamgrv_Parm;  
  // susy parameters  
  std::string m_taudec_Parm; // tauola switch  
  int m_syspin_Parm; // Herwig syspin control  
  int m_spcopt_Parm;  
  int m_itopd; // Top decays switch  

  // Commands to setup herwig from jobOptions  
  CommandVector m_herwigCommandVector;  
  std::vector<long int> m_seeds;

protected:  
  // external proces  int m_ExternalProcess;  
  //*AGB* AcerMC_acset m_acermc_acset;  
  // I/O to HEPEVT  
  //*AGB* void store_Atlas_HEPEVT(void);  
  //*AGB* static Atlas_HEPEVT* atlas_HEPEVT;
};
#endif
