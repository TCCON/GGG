      subroutine compute_tdel(nwin,yobs,yerr,tsen,tdel,terr)
c
c  Computes the temperature error (tdel) and its uncertainty (terr)
c  from column measurements (yobs +/- yerr) made from different
c  windows with different T-sensitivities (tsen).
c
c  Inputs:
c    nwin        I*4  Number of windows
c    yobs(nwin)  R*4  Measured column amounts
c    yerr(nwin)  R*4  Measured column uncertainties
c    tsen(nwin)  R*4  Temperature sensitivity (1/Y.dY/dT)
c
c  Outputs:
c    tdel        R*4  Temperature correction
c    terr        R*4  Temperature uncertainty
c    rchi2       R*4  Reduced CHi2
c    cc          R*4  Correlation Coefficient

c  Subroutine assumes that the retrieved column amounts,
c  Yobs(i) obtained from different windows (same spectrum)
c  have been perturbed from their true value (Ytrue) by a
c  temperature error such that
c      Yobs(i) = Ytrue.[1+Tdel.Tsen(i)]                  (1)
c  where Tsen(i) = 1/Y(i).dY(i)/dT  is the fractional
c  T-sensitivity of the i'th window. There are two
c  unknowns:  Tdel and Tsen.
c
c  Note that this assumes that there are no window-specific
c  biases.  So if the window with the highest E" happens to
c  have a too-low line intensity, the retrieved Tdels
c  will all be biased high.
c
c  This equation can be re-formulated using the same
c  nomenclature as Numerical Recipes:
c      Y(i) = a + b.X(i)                                 (2)
c  where a = Ytrue,  b = Ytrue.Tdel,  and X(i) = Tsen(i)
c  and so Tdel = b/a
c
c  For example if X1=0.5, Y1=25, and X2=0.6, Y2=26,
c    Then       Y=20+10*X
c    or         Y=20*(1+0.5*X)
c  so Tdel = b/a = 10/20 = 0.5
c
c  Subroutine obtains a least-squares estimate of a & b
c  by minimizing
c      Sum_i {[Y(i)-a-b.X(i))]/Yerr(i)}^2                (3)
c  with respect to the two unknowns: a & b. This is
c  achieved my setting their differentials to zero.
c
c     d/da = Sum_i [Y(i)-a-b.X(i))]/Yerr(i)^2 = 0        (4)
c     d/db = Sum_i X(i).[Y(i)-a-b.X(i))]/Yerr(i)^2 = 0   (5)
c
c     Sum Y(i)/Yerr(i)^2 = a.Sum 1/Yerr(i)^2 + b.Sum X(i)/Yerr(i)^2
c     Sum X(i).Y(i)/Yerr(i)^2 = a.Sum X(i)/Yerr(i)^2 + b.Sum X(i)^2/Yerr(i)^2
c
c  resulting in the simultaneous equations
c     Sy  = a.S  + b.Sx                                  (6)
c     Sxy = a.Sx + b.Sxx                                 (7)
c
c  where:
c     S   = Sum_i 1/Yerr(i)^2
c     Sx  = Sum_i X(i)/Yerr(i)^2
c     Sy  = Sum_i Y(i)/Yerr(i)^2
c     Sxy = Sum_i X(i).Y(i)/Yerr(i)^2
c     Sxx = Sum_i X(i)^2/Yerr(i)^2
c
c  D = S.Sxx-Sx^2 is the Denominator or Determinant.
c
c  Multiplying (6) by Sxx and (7) by Sx to eliminate b
c     Sxx.Sy = a.Sxx.S + b.Sxx.Sx
c     Sx.Sxy = a.Sx.Sx + b.Sx.Sxx
c     a = (Sxx.Sy-Sx.Sxy)/(S.Sxx-Sx^2)
c     a = (Sxx.Sy-Sx.Sxy)/D                              (8)
c     ea = SQRT(Sxx/D)                                   (9)
c
c  Multiplying (6) by Sx and (7) by S to eliminate a
c     Sx.Sy  = a.Sx.S  + b.Sx.Sx
c     S.Sxy = a.S.Sx + b.S.Sxx
c     b = (S.Sxy-Sx.Sy)/(S.Sxx-Sx^2)
c     b = (S.Sxy-Sx.Sy)/D                               (10)
c     eb = SQRT(S/D)                                    (11)
c
c  The ea and eb error estimates are independent of
c  the Y-values, they depend only on the x-values and
c  Yerr. They are simply a propagation of the Yerr(i)
c  values through the solution equation and are
c  independent of how well the actual Yobs(i) values
c  conform to a straight line.
c
c  Thus, Tdel = b/a = (S.Sxy-Sx.Sy)/(Sxx.Sy-Sx.Sxy)
c  The fractional uncertainty in Tdel is the RSS of
c  the fractional uncertainties of ea and eb
c  Terr = (b/a).SQRT[(ea/a)^2+(eb/b)^2]
c  Terr = SQRT[((Tdel^2.Sxx+S)/D]/a                     (12)
c  D = S.Sxx-Sx^2 is the Denominator or Determinant.
c
c  In order to prevent zero-divides causing Tdel or Terr= NaN,
c  1.0D-72 is added to the weights (inverse square of Yerr values)
c  to handle situations in which the errors are zero.
c  Also the twxx is set to a minimum value of 1.E-72 to
c  prevent NaNs when the X values (Tsen) are all zero.

      implicit none
      integer*4 nwin,i
      real*4 yobs(nwin),yerr(nwin),tsen(nwin),tdel,terr
      real*8 ww,tw,twx,twy,twxy,twxx,a,denom

      tw  =0.0d0  ! S
      twx =0.0d0  ! Sx
      twy =0.0d0  ! Sy
      twxy=0.0d0  ! Sxy
      twxx=1.0d-72  ! Sxx
      do i=1,nwin
         ww=1.0d0/(1.D-72+dble(yerr(i))**2)
         tw=tw+ww
         twx=twx+ww*tsen(i)
         twy=twy+ww*yobs(i)
         twxy=twxy+ww*tsen(i)*yobs(i)
         twxx=twxx+ww*tsen(i)**2
      end do

c  Thus, Tdel = b/a = (S.Sxy-Sx.Sy)/(Sxx.Sy-Sx.Sxy)
      tdel=sngl((tw*twxy-twx*twy)/(twxx*twy-twx*twxy))

c  Terr = SQRT[((Tdel^2.Sxx+S)/D]/a
      denom=tw*twxx-twx**2
      a=(twxx*twy-twx*twxy)/denom
c      b=(tw*twxy-twx*twy)/denom
      terr=sngl(dsqrt((twxx*tdel**2+tw)/denom)/a)

      return
      end

