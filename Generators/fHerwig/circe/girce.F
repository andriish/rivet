      subroutine girce (x1, x2, p1, p2, rng)
      implicit none
      double precision x1, x2
      integer p1, p2
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
      double precision u, w, circgg
      integer ELECTR, POSITR, PHOTON
      parameter (ELECTR =  11)
      parameter (POSITR = -11)
      parameter (PHOTON =  22)
      if (magic .ne. MAGIC0) then
         call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
      endif
 99   continue
      w = 1d0 / (1d0 + circgg (-1d0, -1d0))
      call rng (u)
      if (u*u .le. w) then
         p1 = POSITR
      else
         p1 = PHOTON
      endif
      call rng (u)
      if (u*u .le. w) then
         p2 = ELECTR
      else
         p2 = PHOTON
      endif
      if (abs(p1) .eq. ELECTR) then
         if (abs(p2) .eq. ELECTR) then
            call gircee (x1, x2, rng)
         elseif (p2 .eq. PHOTON) then
            call girceg (x1, x2, rng)
         endif
      elseif (p1 .eq. PHOTON) then
         if (abs(p2) .eq. ELECTR) then
            call girceg (x2, x1, rng)
         elseif (p2 .eq. PHOTON) then
            call gircgg (x1, x2, rng)
         endif
      endif
      if ((x1 .lt. x1m) .or. (x2 .lt. x2m)) goto 99
      end
