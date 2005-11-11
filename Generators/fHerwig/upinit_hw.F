      SUBROUTINE UPINIT
c     interface to external pythia processes using Les Houches accord
c     Adapted for pythia 6.203
c     Ian Hinchliffe (may 2002)
c     Adapted for Herwig by Borut Paul Kersevan (February 2003)
c 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
 
#include "GeneratorModules/heprup.inc"
#include "GeneratorModules/lesatlas.inc"
      IF(ATLASCHOICE.EQ.3) THEN
         CALL INITACERMC
         CALL ACERMC_ACSET12_HW(-1)
      ENDIF
      IF(ATLASCHOICE.EQ.4) CALL INITALPGEN
      IF(ATLASCHOICE.EQ.5) CALL INITMADGRAPH
      IF(ATLASCHOICE.EQ.6) CALL INITMADCUP
C      IF(ATLASCHOICE.EQ.7) CALL INITTOPREX
      IF(ATLASCHOICE.EQ.8) CALL INITLHAEXT
      IF(ATLASCHOICE.EQ.9) CALL INITMCATNLO
      RETURN
      END
