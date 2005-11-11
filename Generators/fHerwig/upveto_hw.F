      SUBROUTINE UPVETO(IPVETO)
c     interface to external pythia processes using Les Houches accord
c 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
 
#include "GeneratorModules/lesatlas.inc"

      IF(ATLASCHOICE.EQ.4) CALL VETOALPGEN(IPVETO)

      RETURN
      END
