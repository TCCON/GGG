      real*4 function pythag(a,b)
      real*4 a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = max(abs(a),abs(b))
      if (p .eq. 0.0d0) go to 20
      r = (min(abs(a),abs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = sngl(p)
      return
      end

