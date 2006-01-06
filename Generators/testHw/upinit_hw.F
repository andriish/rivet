      SUBROUTINE UPINIT
c     interface to external pythia processes using Les Houches accord
c     Adapted for pythia 6.203
c     Ian Hinchliffe (may 2002)
c     Adapted for Herwig by Borut Paul Kersevan (February 2003)
c 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
 
C#include "GeneratorModules/heprup.inc"
C#include "GeneratorModules/lesatlas.inc"
C      IF(ATLASCHOICE.EQ.3) THEN
C         CALL INITACERMC
C         CALL ACERMC_ACSET12_HW(-1)
C      ENDIF
C      IF(ATLASCHOICE.EQ.4) CALL INITALPGEN
C      IF(ATLASCHOICE.EQ.5) CALL INITMADGRAPH
C      IF(ATLASCHOICE.EQ.6) CALL INITMADCUP
CC      IF(ATLASCHOICE.EQ.7) CALL INITTOPREX
C      IF(ATLASCHOICE.EQ.8) CALL INITLHAEXT
C      IF(ATLASCHOICE.EQ.9) CALL INITMCATNLO
      RETURN
      END
