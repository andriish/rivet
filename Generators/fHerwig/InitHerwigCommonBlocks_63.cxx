//-----------------------------------------------------------------------
//  revision history:
//  -----------------
// *0001 Jan 24 1998 P.Murat: reorganize InitFortranCommonBlocks
// *0002 Oct 29 1998: M. Shapiro:  Change to Herwig 5.9 
// ATLAS VERSION for herwig61/athena Sept 2000 Ian Hinchliffe   
// revised for herwig 6.1 Sept 2001
//-----------------------------------------------------------------------


#include "Herwig_i/herwig63.h"

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
//
Hwumsc_t* gHwumsc;  
				// pointer to Hepevt is defined in 
				// ParticleDB/hepevt.hh
// Hepevt_t* gHepevt;

extern "C" void*  herwig_common_block_address_(char*, int len);

static int Initialized = 0;

					// define pointers to FORTRAN common-blocks
					// it is important that they are returned 
					// by a FORTRAN routine, so it doesn't 
					// require linking in any additional 
					// object files
void InitHerwigCommonBlocks(){

  if (Initialized) return;
  Initialized = 1;
  gHwbeam    = (Hwbeam_t*   ) herwig_common_block_address_("HWBEAM",6); 
  gHwbmch    = (Hwbmch_t*   ) herwig_common_block_address_("HWBMCH",6); 
  gHwproc    = (Hwproc_t*   ) herwig_common_block_address_("HWPROC",6);
  gHwpram    = (Hwpram_t*   ) herwig_common_block_address_("HWPRAM",6);
  gHwprch    = (Hwprch_t*   ) herwig_common_block_address_("HWPRCH",6);
  gHwpart    = (Hwpart_t*   ) herwig_common_block_address_("HWPART",6);
  gHwparp    = (Hwparp_t*   ) herwig_common_block_address_("HWPARP",6);
  gHwbosc    = (Hwbosc_t*   ) herwig_common_block_address_("HWBOSC",6);
  gHwparc    = (Hwparc_t*   ) herwig_common_block_address_("HWPARC",6);
  gHwbrch    = (Hwbrch_t*   ) herwig_common_block_address_("HWBRCH",6);
  gHwevnt    = (Hwevnt_t*   ) herwig_common_block_address_("HWEVNT",6); 
  gHwhard    = (Hwhard_t*   ) herwig_common_block_address_("HWHARD",6); 
  gHwprop    = (Hwprop_t*   ) herwig_common_block_address_("HWPROP",6); 
  gHwunam    = (Hwunam_t*   ) herwig_common_block_address_("HWUNAM",6);  
  gHwupdt    = (Hwupdt_t*   ) herwig_common_block_address_("HWUPDT",6); 
  gHwuwts    = (Hwuwts_t*   ) herwig_common_block_address_("HWUWTS",6); 
  gHwuclu    = (Hwuclu_t*   ) herwig_common_block_address_("HWUCLU",6); 
  gHwdist    = (Hwdist_t*   ) herwig_common_block_address_("HWDIST",6); 
  gHwqdks    = (Hwqdks_t*   ) herwig_common_block_address_("HWQDKS",6); 
  gHwusud    = (Hwusud_t*   ) herwig_common_block_address_("HWUSUD",6);
  // v6.1
  gHwsusy    = (Hwsusy_t*   ) herwig_common_block_address_("HWSUSY",6);
  gHwrpar    = (Hwrpar_t*   ) herwig_common_block_address_("HWRPAR",6);
  gHwminb    = (Hwminb_t*   ) herwig_common_block_address_("HWMINB",6);
  gHwclus    = (Hwclus_t*   ) herwig_common_block_address_("HWCLUS",6);
  // v6.202

  gHw6202    = (Hw6202_t*   ) herwig_common_block_address_("HW6202",6);
  gHwgrav    = (Hwgrav_t*   ) herwig_common_block_address_("HWGRAV",6);
  // v6.301

  gHw6300    = (Hw6300_t*   ) herwig_common_block_address_("HW6300",6);
  gHwcirc    = (Hwcirc_t*   ) herwig_common_block_address_("HWCIRC0",6);
  gHwfmrs    = (Hwfmrs_t*   ) herwig_common_block_address_("HWFMRS",6);
  //
  gHwumsc    = (Hwumsc_t*   ) herwig_common_block_address_("HWUMSC",6); 
  //  gHepevt    = (Hepevt_t*  ) herwig_common_block_address_("HEPEVTD",7); 


}


