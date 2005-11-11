      SUBROUTINE UPEVNT
c     routine for Les Houches interface to external pythir routines
c     adapted from base pythia code. 
c     Ian Hinchliffe (May 2002)
c     Adapted for Herwig by Borut Paul Kersevan (February 2003)

c 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
 
#include "GeneratorModules/hepeup.inc"
#include "GeneratorModules/lesatlas.inc"
      IF(ATLASCHOICE.EQ.3) THEN
         CALL USEACERMC
         CALL ACERMC_ACSET12_HW(0)
      ENDIF
      IF(ATLASCHOICE.EQ.4) CALL USEALPGEN
      IF(ATLASCHOICE.EQ.5) CALL USEMADGRAPH
      IF(ATLASCHOICE.EQ.6) CALL USEMADCUP
C      IF(ATLASCHOICE.EQ.7) CALL USETOPREX
      IF(ATLASCHOICE.EQ.8) CALL USELHAEXT
      IF(ATLASCHOICE.EQ.9) CALL USEMCATNLO
      RETURN
      END
