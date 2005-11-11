      double precision function circeg (x1, x2)
      implicit none
      double precision x1, x2
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
      double precision d1, d2
      if (magic .ne. MAGIC0) then
         call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
      endif
      circeg = -1.0
      if ((ver .eq. 1) .or. (ver .eq. 0)) then
      if (x1 .eq. 1d0) then
         d1 = a1(0)
      elseif (x1 .lt. 1d0 .and. x1 .gt. 0d0) then
         d1 = a1(1) * x1**a1(2) * (1d0 - x1)**a1(3)
      elseif (x1 .eq. -1d0) then
         d1 = 1d0 - a1(0)
      else
         d1 = 0d0
      endif
      if (x2 .lt. 1d0 .and. x2 .gt. 0d0) then
         d2 = a1(4) * x2**a1(5) * (1d0 - x2)**a1(6)
      elseif (x2 .eq. -1d0) then
         d2 = a1(7)
      else
         d2 = 0d0
      endif
      circeg = d1 * d2
      elseif (ver .eq. 2) then
      call circem ('PANIC', '*********************************')
      call circem ('PANIC', '* version 2 has been retired,   *')
      call circem ('PANIC', '* please use version 1 instead! *')
      call circem ('PANIC', '*********************************')
      return
      elseif (ver .gt. 5) then
         call circem ('PANIC', 'versions >5 not available yet')
         return
      else
         call circem ('PANIC', 'version must be positive')
         return
      endif
      end
