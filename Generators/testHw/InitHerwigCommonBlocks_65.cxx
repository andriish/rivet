//-----------------------------------------------------------------------
//  revision history:
//  -----------------
// *0001 Jan 24 1998 P.Murat: reorganize InitFortranCommonBlocks
// *0002 Oct 29 1998: M. Shapiro:  Change to Herwig 5.9 
// ATLAS VERSION for herwig61/athena Sept 2000 Ian Hinchliffe   
// revised for herwig 6.3 Sept 2001
// herwig 6.4 version Jan 2002
//-----------------------------------------------------------------------


#include "Herwig_i/herwig6505.h"

Hwbeam_t* gHwbeam;
Hwbmch_t* gHwbmch;
Hwproc_t* gHwproc;
Hwpram_t* gHwpram;
Hwprch_t* gHwprch;
Hwpart_t* gHwpart;
Hwparp_t* gHwparp;
Hwbosc_t* gHwbosc;
Hwparc_t* gHwparc; 
Hwbrch_t* gHwbrch;
Hwevnt_t* gHwevnt;
Hwhard_t* gHwhard;
Hwprop_t* gHwprop;
Hwunam_t* gHwunam;
Hwupdt_t* gHwupdt;
Hwuwts_t* gHwuwts;
Hwuclu_t* gHwuclu;
Hwdist_t* gHwdist;
Hwqdks_t* gHwqdks;
Hwusud_t* gHwusud;
// 4 new, v61
Hwsusy_t* gHwsusy;
Hwrpar_t* gHwrpar;
Hwminb_t* gHwminb;
Hwclus_t* gHwclus;
// 2 new, v6202
Hwgrav_t* gHwgrav;
Hw6202_t* gHw6202;
// 1 new 6203
Hw6203_t* gHw6203;
// new 6300
Hw6300_t* gHw6300;
Hwcirc_t* gHwcirc;
Hwpmrs_t* gHwpmrs;

// new 6400
Hwdspb_t* gHwdspb;
Hwdsp2_t* gHwdsp2; 
Hwdsp3_t* gHwdsp3;
Hwdsp4_t* gHwdsp4;
Hwdspn_t* gHwdspn;
Hwspin_t* gHwspin;
Hwstau_t* gHwstau;

// new 6500
Hwgupr_t* gHwgupr;
Hw6500_t* gHw6500;
//

// new 6504
Hw6504_t* gHw6504;
//

// Jimmy common blocks
Jmparm_t* gJmparm;

extern "C" void*  herwig_common_block_address__(char*, int len);

static int Initialized = 0;

// Define pointers to FORTRAN common-blocks it is important that they are 
// returned by a FORTRAN routine, so it doesn't require linking in any 
// additional object files
void InitHerwigCommonBlocks() {

  if (Initialized) return;

  gHwbeam    = (Hwbeam_t*   ) herwig_common_block_address__("HWBEAM",6); 
  gHwbmch    = (Hwbmch_t*   ) herwig_common_block_address__("HWBMCH",6); 
  gHwproc    = (Hwproc_t*   ) herwig_common_block_address__("HWPROC",6);
  gHwpram    = (Hwpram_t*   ) herwig_common_block_address__("HWPRAM",6);
  gHwprch    = (Hwprch_t*   ) herwig_common_block_address__("HWPRCH",6);
  gHwpart    = (Hwpart_t*   ) herwig_common_block_address__("HWPART",6);
  gHwparp    = (Hwparp_t*   ) herwig_common_block_address__("HWPARP",6);
  gHwbosc    = (Hwbosc_t*   ) herwig_common_block_address__("HWBOSC",6);
  gHwparc    = (Hwparc_t*   ) herwig_common_block_address__("HWPARC",6);
  gHwbrch    = (Hwbrch_t*   ) herwig_common_block_address__("HWBRCH",6);
  gHwevnt    = (Hwevnt_t*   ) herwig_common_block_address__("HWEVNT",6); 
  gHwhard    = (Hwhard_t*   ) herwig_common_block_address__("HWHARD",6); 
  gHwprop    = (Hwprop_t*   ) herwig_common_block_address__("HWPROP",6); 
  gHwunam    = (Hwunam_t*   ) herwig_common_block_address__("HWUNAM",6);  
  gHwupdt    = (Hwupdt_t*   ) herwig_common_block_address__("HWUPDT",6); 
  gHwuwts    = (Hwuwts_t*   ) herwig_common_block_address__("HWUWTS",6); 
  gHwuclu    = (Hwuclu_t*   ) herwig_common_block_address__("HWUCLU",6); 
  gHwdist    = (Hwdist_t*   ) herwig_common_block_address__("HWDIST",6); 
  gHwqdks    = (Hwqdks_t*   ) herwig_common_block_address__("HWQDKS",6); 
  gHwusud    = (Hwusud_t*   ) herwig_common_block_address__("HWUSUD",6);
  // v6.1
  gHwsusy    = (Hwsusy_t*   ) herwig_common_block_address__("HWSUSY",6);
  gHwrpar    = (Hwrpar_t*   ) herwig_common_block_address__("HWRPAR",6);
  gHwminb    = (Hwminb_t*   ) herwig_common_block_address__("HWMINB",6);
  gHwclus    = (Hwclus_t*   ) herwig_common_block_address__("HWCLUS",6);
  // v6.202

  gHw6202    = (Hw6202_t*   ) herwig_common_block_address__("HW6202",6);
  gHwgrav    = (Hwgrav_t*   ) herwig_common_block_address__("HWGRAV",6);

  // v6.203
  gHw6203    = (Hw6203_t*   ) herwig_common_block_address__("HW6203",6);
  // v6.301

  gHw6300    = (Hw6300_t*   ) herwig_common_block_address__("HW6300",6);
  gHwcirc    = (Hwcirc_t*   ) herwig_common_block_address__("HWCIRC",6);
  gHwpmrs    = (Hwpmrs_t*   ) herwig_common_block_address__("HWPMRS",6);
  // v6.400
  gHwdspb    = (Hwdspb_t*   ) herwig_common_block_address__("HWDSPB",6);
  gHwdsp2    = (Hwdsp2_t*   ) herwig_common_block_address__("HWDSP2",6);
  gHwdsp3    = (Hwdsp3_t*   ) herwig_common_block_address__("HWDSP3",6);
  gHwdsp4    = (Hwdsp4_t*   ) herwig_common_block_address__("HWDSP4",6);
  gHwdspn    = (Hwdspn_t*   ) herwig_common_block_address__("HWDSPN",6);
  gHwspin    = (Hwspin_t*   ) herwig_common_block_address__("HWSPIN",6);
  gHwstau    = (Hwstau_t*   ) herwig_common_block_address__("HWSTAU",6);

  // new 6500
  gHwgupr    = (Hwgupr_t*   ) herwig_common_block_address__("HWGUPR",6);
  gHw6500    = (Hw6500_t*   ) herwig_common_block_address__("HW6500",6);
  // new 6504
  gHw6504    = (Hw6504_t*   ) herwig_common_block_address__("HW6504",6);
  //

  // Jimmy common blocks
  gJmparm    = (Jmparm_t*   ) herwig_common_block_address__("JMPARM",6);;

  Initialized = 1;
}


