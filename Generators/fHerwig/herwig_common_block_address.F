c  For herwig version 6.1 , June-July-2000
      function herwig_common_block_address(common_block_name)
c-----------------------------------------------------------------------
c  returns address of HERWIG common-block by name
c  revision history:
c  -----------------
c *0001 Sep 09 1997 P.Murat: replace 'locf' with 'address' (locf returns 
c                            the address in 4-byte word units)
c  2 Nov 2000 changed to adressi adressc and adressf to remove compile warnings
c   -- Ian Hinchliffe
c-----------------------------------------------------------------------
C      implicit      none
c
c
#include "herwig61.inc"
      character*(*) common_block_name
      integer herwig_common_block_address
c      integer len
c
c-----------------------------------------------------------------------
      integer       adressi
      integer       adressf
      integer       adressc
      integer       common_block_address
c-----------------------------------------------------------------------
      common_block_address = 0
c-----------------------------------------------------------------------
      if (common_block_name.eq.'HWBEAM')then
       herwig_common_block_address = adressi(ipart1)

      else if (common_block_name.eq.'HWBMCH')then
       herwig_common_block_address = adressc(part1)

      else if (common_block_name.eq.'HWPROC')then
       herwig_common_block_address = adressf(ebeam1)

      else if (common_block_name.eq.'HWPRAM')then
       herwig_common_block_address = adressf(afch(1,1))

      else if (common_block_name.eq.'HWPRCH')then
       herwig_common_block_address = adressc(autpdf)

      else if (common_block_name.eq.'HWPART')then
       herwig_common_block_address = adressi(nevpar)
    
      else if (common_block_name.eq.'HWPARP')then
       herwig_common_block_address = adressf(decpar(1,1))

      else if (common_block_name.eq.'HWBOSC')then
       herwig_common_block_address = adressf(alpfac)

      else if (common_block_name.eq.'HWPARC')then
       herwig_common_block_address = adressi(jcopar(1,1))

      else if (common_block_name.eq.'HWBRCH')then
       herwig_common_block_address = adressf(anomsc(1,1))

      else if (common_block_name.eq.'HWEVNT')then
       herwig_common_block_address = adressf(avwgt)
 
      else if (common_block_name.eq.'HWHARD')then
       herwig_common_block_address = adressf(asfixd)

      else if (common_block_name.eq.'HWPROP')then
       herwig_common_block_address = adressf(rltim(0))

      else if (common_block_name.eq.'HWUNAM')then
       herwig_common_block_address = adressc(rname(0))
 
      else if (common_block_name.eq.'HWUPDT')then
       herwig_common_block_address = adressf(brfrac(1))

      else if (common_block_name.eq.'HWUWTS')then
       herwig_common_block_address = adressf(repwt(0,0,0))
 
      else if (common_block_name.eq.'HWUCLU')then
       herwig_common_block_address = adressf(cldkwt(1))
 
      else if (common_block_name.eq.'HWDIST')then
       herwig_common_block_address = adressf(exag)

      else if (common_block_name.eq.'HWQDKS')then
       herwig_common_block_address = adressf(vtxqdk(1,1))

      else if (common_block_name.eq.'HWUSUD')then
       herwig_common_block_address = adressf(accur)     

cc v6.1
      else if (common_block_name.eq.'HWSUSY')then
       herwig_common_block_address = adressf(tanb)
 
      else if (common_block_name.eq.'HWRPAR')then
       herwig_common_block_address = adressf(lamda1(1,1,1))
 
      else if (common_block_name.eq.'HWMINB')then
       herwig_common_block_address = adressf(pmbn1)
 
      else if (common_block_name.eq.'HWCLUS')then
       herwig_common_block_address = adressf(ppcl(1,1))
 
cc
      else if (common_block_name.eq.'HWUMSC')then
       herwig_common_block_address = adressi(lhwout)      

      else  

       write(6,*) " Unknown Common Block Name:  " , common_block_name 

      end if

 
      return
      end

