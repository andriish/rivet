c circe.f -- canonical beam spectra for linear collider physics
c   Copyright (C) 1996,1997 by Thorsten.Ohl@Physik.TH-Darmstadt.de
c
c   Circe is free software; you can redistribute it and/or modify it
c   under the terms of the GNU General Public License as published by
c   the Free Software Foundation; either version 2, or (at your option)
c   any later version.
c
c   Circe is distributed in the hope that it will be useful, but
c   WITHOUT ANY WARRANTY; without even the implied warranty of
c   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c   GNU General Public License for more details.
c
c   You should have received a copy of the GNU General Public License
c   along with this program; if not, write to the Free Software
c   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
      double precision function circe (x1, x2, p1, p2)
      implicit none
      double precision x1, x2
      integer p1, p2
      double precision circee, circeg, circgg
      integer ELECTR, POSITR, PHOTON
      parameter (ELECTR =  11)
      parameter (POSITR = -11)
      parameter (PHOTON =  22)
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
      if (magic .ne. MAGIC0) then
         call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
      endif
      circe = -1.0
      if (abs(p1) .eq. ELECTR) then
         if (abs(p2) .eq. ELECTR) then
            circe = circee (x1, x2)
         elseif (p2 .eq. PHOTON) then
            circe = circeg (x1, x2)
         endif
      elseif (p1 .eq. PHOTON) then
         if (abs(p2) .eq. ELECTR) then
            circe = circeg (x2, x1)
         elseif (p2 .eq. PHOTON) then
            circe = circgg (x1, x2)
         endif
      endif
      end
