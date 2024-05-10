      real*4 function pythag(a,b)
      real*4 a,b
c
c  Finds the length of the hypotenuse of a right-angled triangle
c  whose perpendicular sides are of length a and b.
c     pythag = sqrt(a**2+b**2)
c
c  It does this iteratively without overflow or destructive underflow,
c  without computing any square roots, or squaring a or b. The ratio
c  min(a,b)/max(a,b) is squared, but this can't under-/over-flow.
c
c  Based on:
c  Moler and Morrison, "Replacing Square Roots by Pythagorean Sums", 
c  IBM J. Res. Develop., 27, 1983
c
c  Mathematical Basis.
c  Let  Sqrt(a^2+b^2) = a.Sqrt(1+(b/a)^2) for b<a 
c                     = b.Sqrt(1+(a/b)^2) for a<b 
c                     = p.Sqrt(1+r) 
c
c  Start off with p = max(|a|,|b|)  and  r = (min(|a|,|b|)/p)^2
c
c  The method iteratively reduces r to zero, while increasing p
c  to maintain p.Sqrt(1+r) constant.  When r becomes negligible,
c  p will be the answer we want.
c
c  Approximates Sqrt(1+r) ~ (4+3r)/(4+r) which is accurate to 2nd
c  order for small r. In contrast the classical expansion, 
c  Sqrt(1+r) ~ 1+r/2,  is accurate to first order only.
c
c  Find a new p-value that is close to p.Sqrt(1+r) but does not
c  require the square rooting or squaring. 
c     p' = p.(4+3r)/(4+r)
c
c  Find a new r-value, r', that keeps  p'.Sqrt(1+r') = p.Sqrt(1+r)
c      p.Sqrt(1+r) = p.(4+3r)/(4+r).Sqrt(1+r')
c      1+r' = (1+r).(4+r)^2/(4+3r)^2
c      r' = r^3/(4+3r)^2
c      r' = r.(r/(4+3r))^2
c
c  Then go back to find a new, larger p-value, using r' instead of r.
c
c  In the worst-case scenario of a=b, r=1, and so r'=1/49 = 2.04E-02
c  In the next iteration r' will be reduced to 5.15E-07
c     Iter     p           r
c       0     1.0         1.0     
c       1     1.4         2.04E-002
c       2     1.414213    5.15E-007
c  So the algorithm converges in 3 iterations in the worst-case scenario.

      double precision p,r,s,t,u

      p = max(abs(a),abs(b))
      if (p .gt. 0.0d0) then
         r = (min(abs(a),abs(b))/p)**2
         t = 4.0d0 + r
         do while (t .gt. 4.0d0)
            s = r/t
            u = 1.0d0 + 2.0d0*s
            p = u*p
            r = (s/u)**2 * r
            t = 4.0d0 + r
         end do
      endif
      pythag = sngl(p)
      return
      end
