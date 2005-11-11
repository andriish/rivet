      double precision function girceb (xmin, xmax, a, b, rng)
      implicit none
      double precision xmin, xmax, a, b
      external rng
      double precision t, p, u, umin, umax, x, w
      if ((a .gt. 1d0) .or. (b .lt. 1d0)) then
         girceb = -1d0
         call circem ('ERROR', 'beta-distribution expects a<=1<=b')
         return
      endif
      t = (1d0 - a) / (b + 1d0 - a)
      p = b*t / (b*t + a * (1d0 - t)**b)
      if (xmin .le. 0d0) then
         umin = 0d0
      elseif (xmin .lt. t) then
         umin = p * (xmin/t)**a
      elseif (xmin .eq. t) then
         umin = p
      elseif (xmin .lt. 1d0) then
         umin = 1d0 - (1d0 - p) * ((1d0 - xmin)/(1d0 - t))**b
      else
         umin = 1d0
      endif
      if (xmax .ge. 1d0) then
         umax = 1d0
      elseif (xmax .gt. t) then
         umax = 1d0 - (1d0 - p) * ((1d0 - xmax)/(1d0 - t))**b
      elseif (xmax .eq. t) then
         umax = p
      elseif (xmax .gt. 0d0) then
         umax = p * (xmax/t)**a
      else
         umax = 0d0
      endif
      if (umax .lt. umin) then
         girceb = -1d0
         return
      endif
 10   continue
      call rng (u)
      u = umin + (umax - umin) * u
      if (u .le. p) then
         x = t * (u/p)**(1d0/a)
         w = (1d0 - x)**(b-1d0)
      else
         x = 1d0 - (1d0 - t) * ((1d0 - u)/(1d0 - p))**(1d0/b)
         w = (x/t)**(a-1d0)
      endif
         call rng (u)
      if (w .le. u) goto 10
      girceb = x
      end
