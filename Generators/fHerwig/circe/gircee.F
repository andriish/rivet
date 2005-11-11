      subroutine gircee (x1, x2, rng)
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
      double precision u, girceb
      if (magic .ne. MAGIC0) then
         call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
      endif
      if ((ver .eq. 1) .or. (ver .eq. 0)) then
      call rng (u)
      if (u .le. a1(0)) then
         x1 = 1d0
      else
         x1 = 1d0 - girceb (0d0, 1d0-x1m, a1(3)+1d0, a1(2)+1d0, rng)
      endif
      call rng (u)
      if (u .le. a1(0)) then
         x2 = 1d0
      else
         x2 = 1d0 - girceb (0d0, 1d0-x2m, a1(3)+1d0, a1(2)+1d0, rng)
      endif
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
