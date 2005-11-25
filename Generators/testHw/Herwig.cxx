// ------------------------------------------------------------- 
// File:Herwig.cxx
// Generators/Herwig.cxx Description: Allows the user
// to generate Herwig events and store the result in the
// Transient Store.
//
// AuthorList:
//         Ian Hinchliffe:  Initial Code October : 2000
//         Modeled after the CDF code by Marge Shapriro
//         and Jeremy Lys
//
//
//         Modified for Version 6.4 Jan 2002-- Ian Hinchliffe
// AcerMC 1.x interface added by Borut Paul Kersevan (February 2003)
// TAUOLA/PHOTOS interface added by Borut Paul Kersevan (October 2003)
//
// Header for this module:-

#include "Herwig_i/Herwig.h"
//#include "GeneratorModules/GeneratorName.h"

// Framework Related Headers:-
//#include "GaudiKernel/MsgStream.h"

// Other classes used by this class:-
#include "Herwig_i/herwig6505.h"
#include "HepMC/IO_HERWIG.h"
#include "HepMC/IO_HEPEVT.h"
#include "HepMC/IO_Ascii.h"
#include "HepMC/HEPEVT_Wrapper.h"
//#include "GeneratorUtils/StringParse.h"
//#include <stdlib.h>

#include "CLHEP/Random/RandFlat.h"
//#include "AthenaKernel/IAtRndmGenSvc.h"
//#include "StoreGate/StoreGateSvc.h"

//-------------------------------

#define evtcon evtcon_

extern struct
{
  int istg;
} evtcon;

// Pointer On AtRndmGenSvc
IAtRndmGenSvc* 	p_AtRndmGenSvc;
std::string	herwig_stream	=	"HERWIG_INIT";
extern "C" double atl_hwrgen_( int* idummy )
{
  CLHEP::HepRandomEngine* engine = p_AtRndmGenSvc->GetEngine(herwig_stream);
  return CLHEP::RandFlat::shoot(engine);
}

// calls to fortran routines
extern "C" double atl_phoran_( int* idummy )
{
  CLHEP::HepRandomEngine* engine = p_AtRndmGenSvc->GetEngine("PHOTOS");
  return CLHEP::RandFlat::shoot(engine);
}

extern "C" float atl_ranmar_()
{
  CLHEP::HepRandomEngine* engine = p_AtRndmGenSvc->GetEngine("TAUOLA");
  return (float) CLHEP::RandFlat::shoot(engine);
}

void InitHerwigCommonBlocks();

// calls to fortran routines
extern "C" {
  //  void stdhepctb_();
  // In Herwig_i  area:
  void herwiginterface_(int* itopd);
  //  void cdfreadsusy_(char* file,int *syunit);
  void cdfreadsusy_(const char* ,int * , int);
  //  In Herwig 
  //susy
  void hwissp_(int *syunit);

  // Alpgen special
  //void hwigup_();
    
  void hwuinc_();
  //  void hwusta(int *tstab);
  //  void hwabeg_();
  void hweini_();
  void hwuine_();
  void hwupro_();
  void hwepro_();
  void hwbgen_();
  //void hwdhqk_(); // does not exist in version 6.1
  void hwdhob_();
  void hwcfor_();
  void hwcdec_();
  void hwdhad_();
  void hwdhvy_();
  void hwmevt_();
  void hwufne_(); 
  void hwefin_();
  void hwanal_(int *ihwcod);
  //  void hwefin_();
  //  void hwaend_();
  // In Stdhep
  //    void hwghep_(int *mconv);
  void hwgpdg_(int *mconv);
  void extproc_(int*);

  void upveto_(int*);
  
}
//-----------------------------------
//using HepMC::Vertex;
//using HepMC::Particle;  
using HepMC::IO_HERWIG;
using HepMC::IO_HEPEVT;
using HepMC::IO_Ascii;

// File scope declarations:-

// set pointer to zero at start
Atlas_HEPEVT*  Herwig::atlas_HEPEVT = new Atlas_HEPEVT();

//--------------------------------------------------------------------------
Herwig::Herwig(const std::string& name, 
	       ISvcLocator* pSvcLocator): GenModule(name,pSvcLocator)
{
  //--------------------------------------------------------------------------  
  declareProperty("HerwigCommand", m_herwigCommandVector);
}
//--------------------------------------------------------------------------
Herwig::~Herwig(){
  //--------------------------------------------------------------------------
}
//-------------------------------------------------------------
AcerMC_acset& Herwig::acermc_acset() {
  return m_acermc_acset;
}

//---------------------------------------------------------------------------
StatusCode Herwig::genInitialize() {
  //---------------------------------------------------------------------------
  // Initialise the listing output, parameter and decay data input streams
  //
  MsgStream log(messageService(), name());
  log << MSG:: INFO << " HERWIG INITIALISING.  \n"  << endreq;
  
  static const bool CREATEIFNOTTHERE(true);
  StatusCode RndmStatus = service("AtRndmGenSvc", p_AtRndmGenSvc, CREATEIFNOTTHERE);
  if (!RndmStatus.isSuccess() || 0 == p_AtRndmGenSvc)
    {
      log << MSG::ERROR << " Could not initialize Random Number Service" << endreq;
      return RndmStatus;
    }	

  m_noshower_Parm=false;
  m_nohadroniz_Parm=false;
  // default is to shower and hadronize
  HepMC::HEPEVT_Wrapper::set_sizeof_int(4);
  HepMC::HEPEVT_Wrapper::set_sizeof_real(8);
  HepMC::HEPEVT_Wrapper::set_max_number_entries(10000); 
  // Map Herwig Common blocks to global C structures
  InitHerwigCommonBlocks();

  // Initialize all local parameters to Herwig defaults
  Her2Localpar();  
  
  // set some defaults
  m_Beamtype1_Parm="P       ";
  m_Beamtype2_Parm="P       ";
  m_Energy1_Parm=7000.;  // beam energies
  m_Energy2_Parm=7000.; // beam energy
  m_top_Parm=175.;  // top mass
  m_W_Parm=80.42; // W mass
  m_randomseed1_Parm=17673;
  m_randomseed2_Parm=63565;
  m_process=1400;  // default process is W+Jets
  m_maxer_Parm=500000; // max number of events
  m_modpdf=12;  
  m_ExternalProcess = 0;

  //Set defaults parameters for weak boson decays
  m_modBos_Parm[0] = 0;
  m_modBos_Parm[1] = 0;
  m_modBos_Parm[2] = 0;
  m_modBos_Parm[3] = 0;
  m_modBos_Parm[4] = 0;
  m_modBos_Parm[5] = 0;
  m_modBos_Parm[6] = 0;
  m_modBos_Parm[7] = 0;
  m_modBos_Parm[8] = 0;

  // special loop over the input parameters to load iproc before the call to herwiginterface_()
  CommandVector::const_iterator ic              =       m_herwigCommandVector.begin();
  bool                          proc_not_found  =       true;
  do {
    StringParse mystring(*ic);
    string myvar=mystring.piece(1);
    int myint1=mystring.intpiece(2);
    if(myvar=="iproc") {
      m_process             =       myint1;
      proc_not_found        =       false;
    }
    ic++;
  } while (ic != m_herwigCommandVector.end() && proc_not_found);
  //  initialized  values of Localpars into Herwig Commons
  Localpar2Her();
  // Call Fortran Routine to Initialize Herwig Common Block
  // cout <<  " calling herwiginterface" << endl;
  ic              =       m_herwigCommandVector.begin();
  bool   itopd_not_found  =       true;
  m_itopd = 0;
  do {
    StringParse mystring(*ic);
    string myvar=mystring.piece(1);
    int myint1=mystring.intpiece(2);
    if(myvar=="topdec") {
      m_itopd             =       myint1;
      itopd_not_found        =       false;
    }
    ic++;
  } while (ic != m_herwigCommandVector.end() && itopd_not_found);
  herwiginterface_(&m_itopd);
  // Transfer numbers from Herwig Common Blocks to Localpars
  //   cout <<  " calling Her2Localpar" << endl;
  Her2Localpar();  
  //now read in all the parameter changes from the user
  for (size_t i = 0; i < m_herwigCommandVector.size(); i++) {
    log << MSG:: INFO  << " Command is: " << m_herwigCommandVector[i] << endreq;
    StringParse mystring(m_herwigCommandVector[i]);
    string myvar=mystring.piece(1);
    string myvar1=mystring.piece(2);
    int myint1=mystring.intpiece(2);
    int myint2=mystring.intpiece(3);
    double  myfl1=mystring.numpiece(2);
    double  myfl2=mystring.numpiece(3);
    // Beam types and energy
    if (  myvar== "beam1type") { m_Beamtype1_Parm=myvar1; }
    else if (  myvar== "beam2type") { m_Beamtype2_Parm=myvar1; }
    else if (  myvar== "beam1energy") { m_Energy1_Parm=myfl1; }
    else if (  myvar== "beam2energy") { m_Energy2_Parm=myfl1; }
    // process number and qcd scale
    else if (  myvar== "qcdlam") {m_lambdaQCD_Parm=myfl1;}
    else if(myvar=="iproc"){
      if (myvar1 == "acermc") {
	m_process=-600;
	m_ExternalProcess = 3;
      } else if(myvar1 == "alpgen") {
	m_process=-100;
	m_ExternalProcess = 4;
      } else if (myvar1 == "madgraph") {
	m_process=-950;
	m_ExternalProcess = 5;
      } else if (myvar1 == "madcup") {
	m_process=-800;
	m_ExternalProcess = 6;
      } else if(myvar1 == "lhaext") {
	m_process=-900;
	m_ExternalProcess = 8;
      } else if (myvar1 == "mcatnlo") {
	m_process=-700;
	m_ExternalProcess = 9;
      } else {m_process=myint1;}
      extproc_(&m_ExternalProcess);
    }
    // AcerMC tt~ decay mode switching
    else if(myvar=="acset12") {
      if (m_process==-600) {
	this->acermc_acset().acset12()=myint1;
      }
    }
    // kinematic limits    
    else if(myvar=="ptmin"){m_ptm_Parm=myfl1;}
    else if(myvar=="ptmax"){m_ptM_Parm=myfl1;}
    else if(myvar=="emmin"){m_invarmassm_Parm=myfl1;}
    else if(myvar=="emmax"){m_invarmassM_Parm=myfl1;}
    else if(myvar=="autpdf"){m_autpdf=myvar1;}
    else if(myvar=="modpdf"){m_modpdf=myint1;}
    else if(myvar=="q2min"){ m_q2dilsm_Parm=myfl1;}
    else if(myvar=="pltcut"){ m_pltcut_Parm=myfl1;}
    else if(myvar=="prsof"){ m_prsof_Parm=myfl1;}
    else if(myvar=="clpow"){ m_clpow_Parm=myfl1;}

    // output control and random numbers

    else if(myvar=="iprint"){m_iprint_Parm=myint1;}
    else if(myvar=="maxpr"){m_maxpr_Parm=myint1;}
    else if(myvar=="maxer"){m_maxer_Parm=myint1;}
    else if(myvar=="nrn"){
      if(myint1==1){m_randomseed1_Parm=myint2;}
      else if(myint1==2){ m_randomseed2_Parm=myint2;}
    }
    else if(myvar=="nowgt"){
      if(myvar1=="false"){m_eventweight_Parm=false;}
    }
    else if(myvar=="wgtmax"){m_maxweight_Parm=myfl1;}
    else if(myvar=="effmin"){m_effmin_Parm=myfl1;}
    else if(myvar=="tlout"){m_tlout_Parm=myfl1;}

    // particle masses
    else if(myvar=="rmass"){
      if(myint1==1){m_down_Parm=myfl2;}
      else if(myint1==2){m_up_Parm=myfl2;}
      else if(myint1==3){m_strange_Parm=myfl2;}
      else if(myint1==4){m_charm_Parm=myfl2;}
      else if(myint1==5){m_bottom_Parm=myfl2;}
      else if(myint1==6){m_top_Parm=myfl2;}
      else if(myint1==13){m_gluonmass_Parm=myfl2;}
      else if(myint1==198){m_W_Parm=myfl2;}
      else if(myint1==200){m_Z0_Parm=myfl2;}
      else if(myint1==201){m_higgs_Parm=myfl2;}
      else if(myint1==202){m_ZP_Parm=myfl2;}

    }
    // weak paramters and widths
    else if(myvar=="gamw"){m_wwidth_Parm=myfl1;}
    else if(myvar=="gamz"){m_zwidth_Parm=myfl1;}
    else if(myvar=="gamh"){m_hwidth_Parm=myfl1;}
    else if(myvar=="gamh"){m_hwidth_Parm=myfl1;}
    else if(myvar=="ncolo"){m_colors_Parm=myint1;}
    else if(myvar=="nflav"){m_flavors_Parm=myint1;}
    else if(myvar=="swein"){m_weiangle_Parm=myfl1;}
    else if(myvar=="scabi"){m_cabangle_Parm=myfl1;}

    else if(myvar=="modbos") {
      if(myint1==1){m_modBos_Parm[0]=myint2;}
      else if(myint1==2){m_modBos_Parm[1]=myint2;}
      else if(myint1==3){m_modBos_Parm[2]=myint2;}
      else if(myint1==4){m_modBos_Parm[3]=myint2;}
      else if(myint1==5){m_modBos_Parm[4]=myint2;}
      else if(myint1==6){m_modBos_Parm[5]=myint2;}
      else if(myint1==7){m_modBos_Parm[6]=myint2;}
      else if(myint1==8){m_modBos_Parm[7]=myint2;}
      else if(myint1 >8){m_modBos_Parm[8]=myint2;}
    }

    // susy filename 
    else if(myvar=="susyfile"){m_read_Filesusy=myvar1;}

    // Zprime menu
    else if(myvar=="zprime"){m_zprime_Parm=myint1;}
    else if(myvar=="gamzp"){m_zpwidth_Parm=myfl1;}

    // version 6.3 
    // graviton stuff

    else if(myvar=="grvlam"){m_grvlam_Parm=myfl1;}
    else if(myvar=="emgrv"){m_emgrv_Parm=myfl1;}
    else if(myvar=="gamgrv"){m_gamgrv_Parm=myfl1;}

    // tauola switch
    else if(myvar=="taudec"){m_taudec_Parm=myvar1;}

    // Spin correlations menu
    else if(myvar=="syspin"){m_syspin_Parm=myint1;}
    else if(myvar=="spcopt"){m_spcopt_Parm=myint1;}
    
    else {
      if (myvar!="topdec") log << MSG:: INFO << " ERROR in HERWIG PARAMETERS   " << myvar << " is not a variable name that I recognize. it could be a typo on your part or a mistake on mine. Not all variables are changeable" << endreq;
    }
  }

  // pi0 are NOT Stable
  m_pizstable_Parm = 0;

  // Save the HERWIG_INIT stream seeds....
  CLHEP::HepRandomEngine* engine = p_AtRndmGenSvc->GetEngine(herwig_stream);
  const long*	sip	=	engine->getSeeds();
  long	int	si1	=	sip[0];
  long	int	si2	=	sip[1];

  // Special settings for Alpgen
  if ( m_process == -100 ) m_maxpr_Parm=4;

  //  reset  values of Localpars into Herwig Commons
  Localpar2Her();
  // now for susy
  if(m_process >= 3000 && m_process < 5000){
    int syunit = 66;
    ////    cdfreadsusy_(const_cast<char*>(_read_Filesusy.value().c_str()),&syunit); 
    const string& fileName = m_read_Filesusy; 
    cdfreadsusy_(fileName.c_str(),&syunit,fileName.size()); 
    hwissp_(&syunit);
  }	
  
  // Special initialization for Alpgen
  //  if ( m_process == -100 ) hwigup_();
  
  hwuinc_(); // calculate herwig constants
  hweini_(); // initialise herwig event.
  
  // ... and set them back to the stream for proper save
  p_AtRndmGenSvc->CreateStream(si1, si2, herwig_stream);
  herwig_stream = "HERWIG";

  return StatusCode::SUCCESS;
}

//---------------------------------------------------------------------------
StatusCode Herwig::callGenerator() {
  //---------------------------------------------------------------------------
  MsgStream log(messageService(), name());
  log << MSG:: INFO << " HERWIG generating.  \n"  << endreq;
  const int mxerr = 100;
  int nterr = 0;
  bool goodev = false;

  // Save the random number seeds in the event
  CLHEP::HepRandomEngine*	engine	=	p_AtRndmGenSvc->GetEngine(herwig_stream);
  const long*		s	=	engine->getSeeds();
  m_seeds.clear();
  m_seeds.push_back(s[0]);
  m_seeds.push_back(s[1]);

  while(!goodev){
    hwuine_();  // Initialize event
    if (m_process == -100 ) {
      hwepro_();

      // matching-driven veto
      int ipveto = 0;
      upveto_(&ipveto);
      if(ipveto != 0) {
	gHwevnt->ierror = -1;
	std::cout << " EVENT KILLED BY UPVETO.    EXECUTION CONTINUES" << std::endl;
      }

      if (evtcon.istg <= 0) {
	hwbgen_();  // Generate parton cascade
	hwdhob_();  // Do heavy quark decays 
	hwcfor_();  // Do cluster formation
	hwcdec_();  // Do cluster decays
	hwdhad_();  // Do unstable particle decays
	hwdhvy_();  // Do heavy flavor decays
	hwmevt_();  // Add soft underlying event
      }
    } else {
      hwepro_();  // Generate Hard Process	
      if ( evtcon.istg <= 0 ) {
	if(! m_noshower_Parm) {
	  hwbgen_();  // Generate parton cascade
	  //hwdhqk_();  // Do heavy quark decays removed in version 6.1
	  hwdhob_();  // Do heavy quark decays 
	  if(! m_nohadroniz_Parm) {
	    hwcfor_();  // Do cluster formation
	    hwcdec_();  // Do cluster decays
	    hwdhad_();  // Do unstable particle decays
	    hwdhvy_();  // Do heavy flavor decays
	    hwmevt_();  // Add soft underlying event
	  }
	}
      }
    } 
    hwufne_();  // Finish event
    //  cout << "Event Finished" <<endl;
    // User event analysis if wanted
    if ( evtcon.istg > 0 ) {
      //     if ( HepMC::HEPEVT_Wrapper::number_entries() <= 0 ) {
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      std::cout << "!!! " << generator_name(100000 * m_ExternalProcess) << " TERMINATES NORMALY: NO MORE EVENTS IN FILE !!!" << std::endl;
      std::cout << "!!! PLEASE IGNORE ANY ATHENA ERROR MESSAGES LIKE !!!" << std::endl;
      std::cout << "!!! AthenaEventLoopMgr  ERROR Terminating event processing loop due to errors !!!" << std::endl;
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      return StatusCode::FAILURE;
    }

    int ihwcod = 1;
    hwanal_(&ihwcod);
    if(gHwevnt->ierror > 99 || gHwevnt->ierror < 0) {
      ++nterr;
      if(nterr>mxerr) {
	log << MSG:: INFO << "Herwig:Too many errors on event generation" << endreq;
	goodev=true;
      }
    } else {
      if(ihwcod != 0) {
        goodev=true;
      }
    }
  }

  return StatusCode::SUCCESS;  
  log << MSG:: INFO << " HERWIG generating done.  \n"  << endreq;
}

//---------------------------------------------------------------------------
StatusCode Herwig::genFinalize() {
  //---------------------------------------------------------------------------
  MsgStream log(messageService(), name());
  log << MSG:: INFO << " HERWIG Ending.  \n"  << endreq;
  log << MSG::INFO <<"Call HWEFIN and HWAEND at endRun" << endreq;
  hwefin_();
  // User terminal calculations if wanted
  // hwaend_();
  return StatusCode::SUCCESS;
}
//---------------------------------------------------------------------------
StatusCode Herwig::fillEvt(GenEvent* evt) {
  //---------------------------------------------------------------------------
  MsgStream log(messageService(), name());

  log << MSG::INFO << " HERWIG Atlas_HEPEVT Filing.  \n"  << endreq;
  store_Atlas_HEPEVT();

  log << MSG:: INFO << " HERWIG Filing.  \n"  << endreq;
  HepMC::IO_HERWIG hepio;
  //  hepio.set_trust_mothers_before_daughters(0);
  //hepio.set_print_inconsistency_errors(0);
  //HepMC::HEPEVT_Wrapper::check_hepevt_consistency(); 
  //HepMC::IO_Ascii output("dump.dat",ios::out);
  //  int i = 1;
  // Fill event into HepMC and transient store
  hepio.fill_next_event(evt);
  int pr_id = HERWIG + 100000 * m_ExternalProcess + m_process;
  if (m_process < 0) pr_id = HERWIG + 100000 * m_ExternalProcess;
  if (m_taudec_Parm == "TAUOLA") pr_id += TAUOLA_PHOTOS;
  evt->set_signal_process_id(pr_id);
  evt->set_random_states(m_seeds);
  if (gHwevnt->evwgt < 0.) evt->weights().push_back(-1.);
  if (gHwevnt->evwgt > 0.) evt->weights().push_back(1.);
  if (gHwevnt->evwgt == 0.) log << MSG::WARNING << " EVENT WEIGHT = 0 !!!!!! \n" << endreq;
  //  cout << " ----------------- Print evt " << endl;
  //    evt -> print();
  //  cout << " ----------------- END " << endl;
  
  //  output << evt;
  
  // Convert cm->mm and GeV->MeV
  //   cmTomm(evt);
  GeVToMeV(evt);
  
  return StatusCode::SUCCESS;
}

void Herwig::Her2Localpar() {
  // this changed from CDF code which has .set() and .value() which 
  // are not part of the STL  classes. The belong to the CDF Abspar class.
  string mystring;
  mystring.assign(gHwbmch->part1,8);
  m_Beamtype1_Parm=mystring;
  mystring.assign(gHwbmch->part2,8);
  m_Beamtype2_Parm=mystring;
  mystring.assign(gHwprch->autpdf[0],20);
  m_autpdf=mystring;
  m_modpdf=gHwpram->modpdf[0];
  m_Energy1_Parm=gHwproc ->pbeam1;
  m_Energy2_Parm=gHwproc ->pbeam2;
  m_process=gHwproc->iproc;
  m_maxer_Parm=gHwevnt ->maxer;
  m_randomseed1_Parm=gHwevnt ->nrn[0];
  m_randomseed2_Parm=gHwevnt ->nrn[1];
  m_lambdaQCD_Parm=gHwpram -> qcdlam;
  m_sudord_Parm=gHwusud ->sudord;
  //cjl  _rapiditiesm_Parm=gHwhard ->yjmin;
  m_rapiditiesM_Parm=gHwhard ->yjmax;
  m_ptM_Parm=gHwhard ->ptmax;
  m_ptm_Parm=gHwhard ->ptmin;
  m_invarmassM_Parm=gHwhard ->emmax;
  m_invarmassm_Parm=gHwhard ->emmin;
  m_thrust_Parm=gHwhard ->thmax;
  m_ptpow_Parm=gHwhard ->ptpow;
  m_empow_Parm=gHwhard ->empow;
  m_q2pow_Parm=gHwhard ->q2pow;
  m_qlim_Parm=gHwhard ->qlim;
  m_q2dilsM_Parm=gHwhard ->q2max;
  m_q2dilsm_Parm=gHwhard ->q2min;
  m_bgshat_Parm=gHwhard ->bgshat;

  m_pltcut_Parm=gHwdist ->pltcut;
  
  m_azsoft_Parm=gHwpram->azsoft;
  m_prsof_Parm=gHwpram->prsof;
  m_azspin_Parm=gHwpram -> azspin;
  m_qspac_Parm=gHwpram -> qspac;
  m_ispac_Parm=gHwpram -> ispac;
  m_nospac_Parm=gHwpram -> nospac;
  m_vqcut_Parm=gHwpram -> vqcut;
  m_vgcut_Parm=gHwpram -> vgcut;
  m_vpcut_Parm=gHwpram -> vpcut;
  m_intrinsicpt_Parm=gHwpram -> ptrms;

  m_clmax_Parm=gHwpram -> clmax; 
  m_clpow_Parm=gHwpram -> clpow;
  m_psplt1_Parm=gHwpram -> psplt[0];
  m_psplt2_Parm=gHwpram -> psplt[1];
  m_cldir1_Parm=gHwpram -> cldir[0];
  m_cldir2_Parm=gHwpram -> cldir[1];
  m_clsmr1_Parm=gHwpram -> clsmr[0];
  m_clsmr2_Parm=gHwpram -> clsmr[1];
  m_qdiqk_Parm=gHwpram -> qdiqk;
  m_pdiqk_Parm=gHwpram -> pdiqk;
  m_decwt_Parm=gHwuwts ->decwt;
  m_pwt_Parm=gHwuwts ->pwt[0];

  //cjl need - sign on pmbn3 and pmbk2
  m_pmbn1_Parm=gHwminb -> pmbn1;
  m_pmbn2_Parm=gHwminb -> pmbn2;
  m_pmbn3_Parm=-gHwminb -> pmbn3;
  m_pmbk1_Parm=gHwminb -> pmbk1;
  m_pmbk2_Parm=-gHwminb -> pmbk2;
  m_pmbm1_Parm=gHwminb -> pmbm1;
  m_pmbm2_Parm=gHwminb -> pmbm2;
  m_pmbp1_Parm=gHwminb -> pmbp1;
  m_pmbp2_Parm=gHwminb -> pmbp2;
  m_pmbp3_Parm=gHwminb -> pmbp3;
  m_ensof_Parm=gHwpram -> ensof;

  m_pizstable_Parm=gHwupdt -> rstab[21];

  m_eventweight_Parm=gHwevnt ->nowgt;
  m_maxweight_Parm=gHwevnt ->wgtmax;
  m_effmin_Parm=gHwpram ->effmin;
  m_tlout_Parm=gHwevnt ->tlout;

  m_colors_Parm=gHwpram -> ncolo;
  m_flavors_Parm=gHwpram -> nflav; 
  m_weiangle_Parm=gHwpram -> swein;
  m_cabangle_Parm=gHwpram -> scabi;
  m_wwidth_Parm=gHwpram -> gamw;
  m_zwidth_Parm=gHwpram -> gamz;
  m_hwidth_Parm=gHwpram -> gamh;

  m_down_Parm=gHwprop ->rmass[1];
  m_up_Parm=gHwprop ->rmass[2];
  m_strange_Parm=gHwprop ->rmass[3];
  m_charm_Parm=gHwprop ->rmass[4];
  m_bottom_Parm=gHwprop ->rmass[5]; 
  m_top_Parm=gHwprop ->rmass[6];
  m_gluonmass_Parm=gHwprop ->rmass[13];
  m_W_Parm=gHwprop ->rmass[198];
  m_Z0_Parm=gHwprop ->rmass[200];
  m_higgs_Parm=gHwprop ->rmass[201];
 
  m_modBos_Parm[0] = gHwbosc->modbos[0];
  m_modBos_Parm[1] = gHwbosc->modbos[1];
  m_modBos_Parm[2] = gHwbosc->modbos[2];
  m_modBos_Parm[3] = gHwbosc->modbos[3];
  m_modBos_Parm[4] = gHwbosc->modbos[4];
  m_modBos_Parm[5] = gHwbosc->modbos[5]; 
  m_modBos_Parm[6] = gHwbosc->modbos[6]; 
  m_modBos_Parm[7] = gHwbosc->modbos[7]; 
  m_modBos_Parm[8] = gHwbosc->modbos[8];
 
  m_zprime_Parm=gHwpram -> zprime;
  m_ZP_Parm=gHwprop -> rmass[202];
  m_zpwidth_Parm=gHwpram -> gamzp;

  m_iprint_Parm=gHwpram -> iprint;
  m_maxpr_Parm=gHwevnt -> maxpr;
  m_prvtx_Parm=gHwpram -> prvtx;
  m_prndef_Parm=gHwpram -> prndef;
  m_prntex_Parm=gHwpram -> prntex;
  m_prnweb_Parm=gHwpram -> prnweb;
  m_nprfmt_Parm=gHwpram -> nprfmt;

  m_grvlam_Parm=gHwgrav ->grvlam;
  m_emgrv_Parm=gHwgrav ->emgrv;
  m_gamgrv_Parm=gHwgrav ->gamgrv;

  mystring.assign(gHwdspn->taudec,6);
  m_taudec_Parm=mystring;

  m_syspin_Parm=gHwdspn -> syspin;
  m_spcopt_Parm=gHwspin -> spcopt;


  // Following two are used only internally to this routine:
  //_noshower_Parm=0);// correspondence in common block not finded
  //_nohadroniz_Parm=0);// correspondence in common block not finded
}

//
// This routine takes the Localpars of the Herwig  and transfers them
// to the Herwig Common Blocks

void Herwig::Localpar2Her() {
  int nch;
  nch = (m_Beamtype1_Parm.size()<8) ? m_Beamtype1_Parm.size() : 7;
  for(int i=0; i<8; i++) {
    if(i<nch) {
      gHwbmch->part1[i] = m_Beamtype1_Parm.data()[i];
    }
    else {
      gHwbmch->part1[i] = ' ';
    }
  }
  nch = (m_Beamtype2_Parm.size()<8) ? m_Beamtype2_Parm.size() : 7;
  for(int i=0; i<8; i++) {
    if(i<nch) {
      gHwbmch->part2[i] = m_Beamtype2_Parm.data()[i];
    }
    else {
      gHwbmch->part2[i] = ' ';
    }
  }
  nch = (m_autpdf.size()<20) ? m_autpdf.size() : 19;
  for(int i=0; i<20; i++) {
    if(i<nch) {
      gHwprch->autpdf[0][i] = m_autpdf.data()[i];
      gHwprch->autpdf[1][i] = m_autpdf.data()[i];
    }
    else {
      gHwprch->autpdf[0][i] = ' ';
      gHwprch->autpdf[1][i] = ' ';
    }
  }

  gHwpram ->modpdf[0]=m_modpdf;
  gHwpram ->modpdf[1]=m_modpdf;
  gHwproc ->pbeam1=m_Energy1_Parm;
  gHwproc ->pbeam2=m_Energy2_Parm;
  gHwproc ->iproc=m_process;
  gHwevnt ->maxer=m_maxer_Parm;
  gHwevnt ->nrn[0]=m_randomseed1_Parm;
  gHwevnt ->nrn[1]=m_randomseed2_Parm;
  gHwpram -> qcdlam=m_lambdaQCD_Parm;
  gHwusud ->sudord=m_sudord_Parm;

  //cjl re ymin
  //gHwhard ->yjmin=_rapiditiesm_Parm;
  float yjmin1;
  yjmin1 = -m_rapiditiesM_Parm;
  //cout << " yjmin in Localpar2Her in HM.cc  " << yjmin1  <<endl;
  gHwhard ->yjmin=yjmin1;
  //cjl_e
  gHwhard ->yjmax=m_rapiditiesM_Parm;
  gHwhard ->ptmin=m_ptm_Parm;
  gHwhard ->ptmax=m_ptM_Parm;
  gHwhard ->emmin=m_invarmassm_Parm;
  gHwhard ->emmax=m_invarmassM_Parm;
  gHwhard ->thmax=m_thrust_Parm;
  gHwhard ->ptpow=m_ptpow_Parm;
  gHwhard ->empow=m_empow_Parm;
  gHwhard ->q2pow=m_q2pow_Parm;
  gHwhard ->qlim=m_qlim_Parm;
  gHwhard ->q2min=m_q2dilsm_Parm;
  gHwhard ->q2max=m_q2dilsM_Parm;
  gHwhard ->bgshat=m_bgshat_Parm;

  gHwdist ->pltcut=m_pltcut_Parm;

  gHwpram ->azsoft=m_azsoft_Parm;
  gHwpram ->prsof=m_prsof_Parm;
  gHwpram -> azspin=m_azspin_Parm;
  gHwpram -> qspac=m_qspac_Parm;
  gHwpram -> ispac=m_ispac_Parm;
  gHwpram -> nospac=m_nospac_Parm;
  gHwpram -> vqcut=m_vqcut_Parm;
  gHwpram -> vgcut=m_vgcut_Parm;
  gHwpram -> vpcut=m_vpcut_Parm;
  gHwpram -> ptrms=m_intrinsicpt_Parm;

  gHwpram -> clmax=m_clmax_Parm; 
  gHwpram -> clpow=m_clpow_Parm;
  gHwpram -> psplt[0]=m_psplt1_Parm;
  gHwpram -> psplt[1]=m_psplt2_Parm;
  gHwpram -> cldir[0]=m_cldir1_Parm;
  gHwpram -> cldir[1]=m_cldir2_Parm;
  gHwpram -> clsmr[0]=m_clsmr1_Parm;
  gHwpram -> clsmr[1]=m_clsmr2_Parm;
  gHwpram -> qdiqk=m_qdiqk_Parm;
  gHwpram -> pdiqk=m_pdiqk_Parm;
  gHwuwts -> decwt=m_decwt_Parm;
  gHwuwts -> pwt[0]=m_pwt_Parm;

  // cjl note pmbn3 and pmbk2 negatives here
  gHwminb -> pmbn1=m_pmbn1_Parm;
  gHwminb -> pmbn2=m_pmbn2_Parm;
  gHwminb -> pmbn3=-m_pmbn3_Parm;
  gHwminb -> pmbk1=m_pmbk1_Parm;
  gHwminb -> pmbk2=-m_pmbk2_Parm;
  gHwminb -> pmbm1=m_pmbm1_Parm;
  gHwminb -> pmbm2=m_pmbm2_Parm;
  gHwminb -> pmbp1=m_pmbp1_Parm;
  gHwminb -> pmbp2=m_pmbp2_Parm;
  gHwminb -> pmbp3=m_pmbp3_Parm;
  gHwpram -> ensof=m_ensof_Parm;

  gHwupdt -> rstab[21]=m_pizstable_Parm;

  gHwevnt -> nowgt=m_eventweight_Parm;
  gHwevnt -> wgtmax=m_maxweight_Parm;
  gHwpram -> effmin=m_effmin_Parm;
  gHwevnt -> tlout=m_tlout_Parm;
  gHwpram -> ncolo=m_colors_Parm;
  gHwpram -> nflav=m_flavors_Parm; 
  gHwpram -> swein=m_weiangle_Parm;
  gHwpram -> scabi=m_cabangle_Parm;
  gHwpram -> gamw=m_wwidth_Parm;
  gHwpram -> gamz=m_zwidth_Parm;
  gHwpram -> gamh=m_hwidth_Parm;

  gHwprop ->rmass[1]=m_down_Parm;
  gHwprop ->rmass[2]=m_up_Parm;
  gHwprop ->rmass[3]=m_strange_Parm;
  gHwprop ->rmass[4]=m_charm_Parm;
  gHwprop ->rmass[5]=m_bottom_Parm; 
  gHwprop ->rmass[6]=m_top_Parm;
  gHwprop ->rmass[13]=m_gluonmass_Parm;
  gHwprop ->rmass[198]=m_W_Parm;
  gHwprop ->rmass[200]=m_Z0_Parm;
  gHwprop ->rmass[201]=m_higgs_Parm;

  for ( int i = 0; i<9; i++ )
    gHwbosc->modbos[i] = m_modBos_Parm[i];

  gHwpram -> zprime = m_zprime_Parm;
  gHwprop -> rmass[202] = m_ZP_Parm;
  gHwpram -> gamzp = m_zpwidth_Parm;

  gHwpram -> iprint = m_iprint_Parm;
  gHwevnt -> maxpr = m_maxpr_Parm;
  gHwpram -> prvtx = m_prvtx_Parm;
  gHwpram -> prndef = m_prndef_Parm;
  gHwpram -> prntex = m_prntex_Parm;
  gHwpram -> prnweb = m_prnweb_Parm;
  gHwpram -> nprfmt = m_nprfmt_Parm;

  gHwgrav -> grvlam = m_grvlam_Parm;
  gHwgrav -> emgrv = m_emgrv_Parm;
  gHwgrav -> gamgrv = m_gamgrv_Parm;

  nch = (m_taudec_Parm.size()<7) ? m_taudec_Parm.size() : 6;
  for(int i=0; i<6; i++) {
    if(i<nch) {
      gHwdspn->taudec[i] = m_taudec_Parm.data()[i];
    }
  }

  gHwdspn -> syspin = m_syspin_Parm;
  gHwspin -> spcopt = m_spcopt_Parm;

}

void
Herwig::store_Atlas_HEPEVT(void)
{
  MsgStream log(messageService(), name());

  // std::cout << "atlas_HEPEVT------" << atlas_HEPEVT->nhep()  << std::endl;
  //std::cout << "atlas_HEPEVT------" << atlas_HEPEVT->isthep(10)  << std::endl;
  //std::cout << "atlas_HEPEVT------" << atlas_HEPEVT->idhep(10)  << std::endl;
  //std::cout << "atlas_HEPEVT------" << atlas_HEPEVT->jmohep(1,10)  << std::endl;
  //std::cout << "atlas_HEPEVT------" << atlas_HEPEVT->jdahep(2,10)  << std::endl;

  atlas_HEPEVT->fill();

  Atlas_HEPEVT* Ahep = new Atlas_HEPEVT();
  *(Ahep)=*(atlas_HEPEVT);
  std::string keyid = "Herwig";
  m_sgSvc->record(Ahep, keyid);

}
