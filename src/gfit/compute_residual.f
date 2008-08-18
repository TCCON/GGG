      subroutine compute_residual(mode,obsrvd,calcul,resids,nmp)
c
c  Computes the difference between vectors OBSRVD and CALCUL
c  If MODE=1 it returns the arithmetic residual = OBSRVD-CALCUL
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
c  The problem with so-called logarithmic residuals is that
c  the LOG function has a nasty habit of returning NAN's if,
c  for example, OBSRVD and CALCUL are opposite in sign or if
c  either of them is zero. This can happen due to measurement
c  noise, or if a spectral region is strongly blacked out.
c
c  By contrast, the usual method of defining residuals as
c  OBSRVD-CALCUL has no possible problem with NaN's
c
c  So how can we take advantage of the fast convergence
c  afforded by the logarithmic residuals, without becoming
c  hamstrung by the occasional NaN's?
c
c  Define del=(OBSRVD-CALCUL)/CALCUL
c  For good spectral fits OBSRVD ~ CALCUL and so del is small
c  Exploiting the approximation LOGe(1+del) ~ del for small del
c    RESIDUAL=CALCUL*LOGe(OBSRVD/CALCUL)
c            =CALCUL*LOGe(1+(OBSRVD-CALCUL)/CALCUL)
c            =CALCUL*(OBSRVD-CALCUL)/CALCUL)
c            =OBSRVD-CALCUL
c
c This means that if the spectral fits are good, it makes no
c difference whether we use arithmetic or logarithmic residuals 
c They are the same thing.
c
c  To avoid this problem we don't actually use the LOGe function.
c  Instead we approximate LOGe(1+del)= del - del^2/2 + del^3/3
c  LOGe(1+del)= del*(1-0.5*del*(1+0.6666*del) 
c  CALCUL*LOGe(1+del)= RES*(1-0.5*(RES/CALCUL)*(1+0.666*RES/CALCUL)) 
c    RESIDUAL=CALCUL*LOGe(OBSRVD/CALCUL)
c            =CALCUL*LOGe(1+(OBSRVD-CALCUL)/CALCUL)
c            =CALCUL*(OBSRVD-CALCUL)/CALCUL*(1-0.5*(OBSRVD-CALCUL)/CALCUL)
c            =(OBSRVD-CALCUL)*(1-0.5*(OBSRVD-CALCUL)/CALCUL)
c            =(OBSRVD-CALCUL)*(1-0.5*(OBSRVD/CALCUL-1))
c            =(OBSRVD-CALCUL)*(1-0.5*OBSRVD/CALCUL+0.5)
c            =(OBSRVD-CALCUL)*(1.5-0.5*OBSRVD/CALCUL)
c            =(OBSRVD-CALCUL)*(3-OBSRVD/CALCUL)/2
c            =(OBSRVD-CALCUL)*(2*CALCUL+CALCUL-OBSRVD)/CALCUL/2
c   RES= OBSRVD-CALCUL
c   CORR = RES/CALCUL
c   CALCUL*LOGe(1+del)= RES*(1-0.5*CORR*(1+0.666*CORR)) 

      integer*4 mode,nmp,imp
      real*4  obsrvd(nmp),calcul(nmp),resids(nmp)

c  Calculate residuals
      call vsub(obsrvd,1,calcul,1,resids,1,nmp) ! residuals
      if(mode.gt.1) then ! use logarithmic residuals 1'st convergence
         do imp=1,nmp
            resoc=resids(imp)/calcul(imp)
            if(abs(resoc).lt.0.5) resids(imp)=resids(imp)*(1-resoc/2)
         end do   ! imp=1,nmp
      endif   !  mode.gt.1
      return
      end
