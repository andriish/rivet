      subroutine circem (errlvl, errmsg)
      implicit none
      character*(*) errlvl, errmsg
      integer MAGIC0
      parameter (MAGIC0 = 1904 06 16)
      double precision x1m, x2m, roots
      common /circom/ x1m, x2m, roots
      double precision lumi
      common /circom/ lumi
      double precision a1(0:7)
      common /circom/ a1
      integer acc, ver, rev, chat
      common /circom/ acc, ver, rev, chat
      integer magic
      common /circom/ magic
      integer e, r
      common /circom/ e, r
      save /circom/
      integer errcnt
      save errcnt
      data errcnt /0/
      if (errlvl .eq. 'MESSAGE') then
         print *, 'circe:message: ', errmsg
      elseif (errlvl .eq. 'WARNING') then
         if (errcnt .lt. 100) then
            errcnt = errcnt + 1
            print *, 'circe:warning: ', errmsg
         elseif (errcnt .eq. 100) then
            errcnt = errcnt + 1
            print *, 'circe:message: more than 100 messages'
            print *, 'circe:message: turning warnings off'
         endif
      elseif (errlvl .eq. 'ERROR') then
         if (errcnt .lt. 200) then
            errcnt = errcnt + 1
            print *, 'circe:error:   ', errmsg
         elseif (errcnt .eq. 200) then
            errcnt = errcnt + 1
            print *, 'circe:message: more than 200 messages'
            print *, 'circe:message: turning error messages off'
         endif
      elseif (errlvl .eq. 'PANIC') then
         if (errcnt .lt. 300) then
            errcnt = errcnt + 1
            print *, 'circe:panic:   ', errmsg
         elseif (errcnt .eq. 300) then
            errcnt = errcnt + 1
            print *, 'circe:message: more than 300 messages'
            print *, 'circe:message: turning panic messages off'
         endif
      else
         print *, 'circe:panic:    invalid error code ', errlvl
      endif
      end
