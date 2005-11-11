      subroutine gircgg (x1, x2, rng)
      implicit none
      double precision x1, x2
      external rng
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
      double precision girceb
      if (magic .ne. MAGIC0) then
         call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
      endif
      if ((ver .eq. 1) .or. (ver .eq. 0)) then
      x1 = girceb (x1m, 1d0, a1(5)+1d0, a1(6)+1d0, rng)
      x2 = girceb (x2m, 1d0, a1(5)+1d0, a1(6)+1d0, rng)
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
