      PROGRAM fTest
        CALL writeIt()
        CALL writeThat(12,32.67896d0,'d')
      END PROGRAM fTest

      SUBROUTINE writeIt()
        implicit none
        character*20 foo
        foo = "hello world"
        WRITE(*,*) foo
      END SUBROUTINE writeIt

      SUBROUTINE writeThat(bar,baz,whee)
        implicit none
        integer bar
        double precision baz
        character whee
10      FORMAT (I3,X,D10.5,X,A1)
        WRITE(*,10) bar, baz, whee
      END SUBROUTINE writeThat
