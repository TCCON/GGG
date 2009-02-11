      subroutine compute_residual(mode,obsrvd,calcul,resids,nmp)
c
c  Computes the difference between vectors OBSRVD and CALCUL
c  If MODE=1 it returns the arithmetic residual = OBSRVD - CALCUL
c  If MODE=2 it returns the logarithmic residual = CALCUL*LOGe[OBSRVD/CALCUL]
c
c  Inputs:
c        MODE         I*4   Mode
c        OBSRVD(NMP)  R*4   The observed/measured spectrum
c        CALCUL(NMP)  R*4   The forward_model/calculated spectrum
c        NMP          I*4   Number of Measured Points
c
c  Output:
c        RESIDS(NMP)  R*4   The residual spectrum
c
c  Explanation:
c     In situations where you are solving a non-linear equation
c  of the form  y(i) = a.exp(-c(i).x(j)), you achieve much
c  faster convergence if you calculate the residuals as
c     RESIDS=CALCUL*LOGe(OBSRVD/CALCUL)
c  rather than the more usual
c     RESIDS=OBSRVD-CALCUL
c
c  To understand why, assume that
c  x' is the current best estimate of the state vector, and
c  x is its unknown true value. Hence
c     OBSRVD = a.exp(-c.x)
c     CALCUL = a.exp(-c.x')
c  Therefore
c     OBSRVD/CALCUL = exp(-c.(x-x'))
c     LOGe[OBSRVD/CALCUL] = -c.(x-x')
c  So   x = x' - LOGe[OBSRVD/CALCUL]/c 
c  So the update to x' is proportional to LOGe[OBSRVD/CALCUL]
c
c  The problem with so-called logarithmic residuals is that the
c  LOGe function has a nasty habit of returning NAN's or Inf's
c  if, for example, OBSRVD and CALCUL are opposite in sign or if
c  either of them is zero. This can happen due to measurement
c  noise, or if a spectral region is strongly blacked out.
c
c  By contrast, the usual method of defining residuals as
c  OBSRVD-CALCUL has no possible problem with NaN's or Inf's.
c
c  So how can we take advantage of the fast convergence
c  afforded by the logarithmic residuals, without becoming
c  hamstrung by the occasional NaN's?
c
c  To avoid this problem we don't actually use the LOGe function.
c  We use a power series expansion of the LOGe function , i.e.
c     LOGe(1+x) = x - 0.5*x^2 
c
c  Define RES = OBSRVD-CALCUL
c
c    LOGe[OBSRVD/CALCUL] = LOGe[1+RES/CALCUL]
c                        = RES/CALCUL - 0.5*(RES/CALCUL)**2 
c
c  So CALCUL*LOGe[OBSRVD/CALCUL] = RES*[1-0.5*(RES/CALCUL)]
c
c  This means that if the spectral fits are good (RES/CALCUL << 1)
c  it makes no difference whether we use arithmetic or logarithmic
c  residuals -- they are the same thing.
c
      integer*4 mode,nmp,imp
      real*4  obsrvd(nmp),calcul(nmp),resids(nmp)

c  Calculate residuals
      call vsub(obsrvd,1,calcul,1,resids,1,nmp) ! residuals
      if(mode.gt.1) then ! use logarithmic residuals 1'st convergence
         do imp=1,nmp
c            resids(imp)=obsrvd(imp)-calcul(imp)
            resoc=resids(imp)/calcul(imp)
            if(abs(resoc).lt.0.5) resids(imp)=resids(imp)*(1-resoc/2)
         end do   ! imp=1,nmp
      endif   !  mode.gt.1
      return
      end
