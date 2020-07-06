      subroutine siv_correction(nip, margin, bufigram, bufsmoo,
     &  frsp, pinl, pinv, sivcflag, dclevel, fvsi_calc, frzpda)

c  Corrects an interferogram for Solar Intensity Variations (SIV) by dividing
c  the raw DC interferogram (BUFIGRAM) by a smoothed version of itself.
c  This restores fringes affected by SIVs to their true amplitudes,
c  thereby producing a corrected interferogram that represents what
c  would have been measured under conditions of constant luminosity. 
c
c  The subroutine also determines the peak interferogram location (PINL)
c  using the uncorrected AC interferogram values.
c 
c  Inputs:
c           nip           I*4  Number of interferogram points
c           margin        I*4  Number of points at each end of igram
c           bufigram(nip) R*4  Raw Interferogram (AC or DC)
c           bufsmoo(nip)  R*4  Workspace for smoothed igram
c           frsp          R*4  Fraction of spectrum passed by filter
c
c Outputs: 
c           bufigram(nip) R*4  Corrected AC Interferogram
c           pinl          I*4  ZPD location
c           pinv          I*4  Interferogram  Value at ZPD
c           sivcflag      I*4  SIV Correction flag
c           dclevel       R*4  DC-level at ZPD
c           frzpda        R*4  Fractional Depth of dip at ZPD 
c                              (DIP in spectrum header)
c
c
c  If the input igram is AC, the subroutine will simply locate ZPD and
c  subtract the mean igram level.
c
c  If the input igram is DC, but an frsp-value <= 0 is specified in the
c  i2s.in file, the subroutine will subtract the smoothed igram from
c  the raw igram, but not divide the result by the smoothed igram
c
c  On output, BUFIGRAM always contains an AC igram, irrespective of
c  whether the raw igram was AC or DC. This simplifies the subsequent
c  code (e.g. i2s_lite_siv.f).
c
c  Smoothing is performed by convolution with a smoothing operator
c  designed to pass no signal at frequencies > FRSP*Nyquist.
c  Smoothed interferogram is modified around ZPD to remove any 
c  artifacts (e.g. dip associated with detection non-linearity).
c
c  This subroutine has four main differences from Gretchen Keppel-Alek's
c  original FFT-based dc-correction subroutine (Applied Optics, 2007): 
c  1) More robust in terms of finding the true ZPD
c  2) Conceptual simplicity: we're simply convolving the
c     interferogram with a smoothing operator. No FFTs.
c  3) It is faster for typical operator lengths < 2.log2(NFFT)
c  4) The smoothing is performed in a separate subroutine, 
c     which also fits a straight line is fitted through the
c     ZPD region of the smoothed igram, avoiding any distortion
c     due to detection non-linearity or slew-rate limitations.
c  
c  This subroutine has an improved method for finding ZPD.
c  Instead of using the uncorrected DC igram values or the
c  corrected DC igram values, it uses the uncorrected AC igram values.
c  The former methods can produce false ZPD's when the
c  DC interferogram falls close to zero (e.g. clouds).
c  Such igrams are not necessarily useless, for example,
c  the low-signal occuring on the short side of ZPD.
c
      implicit none
      integer*4
     & margin,kmin,kmax,sivcflag,
     & nip,jip,pinl,nglitch

      real*4
     & bufigram(nip), ! Data storage and FFT bufigram
     & bufsmoo(nip), ! Smoothed inteferogram
     & ymin,ymax,ybar,
     & smoo_min,smoo_max,soffset,sthresh,smoo_mean,smoo_std,
     & pinv,frsp,gthresh

      parameter (sthresh=0.01,gthresh=1.0E+18)

      real*8 dclevel,fvsi_calc,frzpda,tot

c  Determine min, max, and mean of raw interferogram values,
c  skipping ends of igram that won't be phase corrected.
c  Replace single-point glitches by the mean of its neighbours.
      kmin=1+margin
      kmax=1+margin
      ymax=bufigram(kmin)
      ymin=bufigram(kmax)
      ybar=bufigram(kmin)
      nglitch=0
      do jip=1+margin,nip-margin
         if(abs(bufigram(jip)).gt.gthresh) then
            bufigram(jip)=0.5*(bufigram(jip-1)+bufigram(jip+1))
            if(abs(bufigram(jip)).gt.gthresh) then
               bufigram(jip)=sign(gthresh,bufigram(jip))
            endif
            nglitch=nglitch+1
         endif
         ybar=ybar+bufigram(jip)
         if(bufigram(jip).lt.ymin) then
            kmin=jip
            ymin=bufigram(jip)
         endif
         if(bufigram(jip).gt.ymax) then
            kmax=jip
            ymax=bufigram(jip)
         endif
      end do
      ybar=ybar/nip
      if(nglitch.ge.1)write(*,*)'Warning! Corrected ',
     & nglitch,' glitch(es).'
c      write(*,*) 'siv_correction: kmin,kmax=',kmin,kmax
c      write(*,*) 'siv_correction: ymin,ymax=',ymin,ymax
c      write(*,*) 'siv_correction: ybar=',ybar
 
c  Determine whether interferogram is AC- or DC-coupled.
c  This is done from the interferogram itself, not the header info.
c  If all the igram values have the same sign (ymin*ymax>0) it is DC.
c  But if ymin and ymax have opposite sign, igram could be AC or DC.
c  So compare Ymax-Ymin with Ybar.
c  If AC, determine PINL, PINV, subtract mean value and return
      if( abs(ybar).lt.0.2*(ymax-ymin) ) then      ! Treat as AC
         write(*,*)'AC Interferogram: ymin,ybar,ymax ', ymin,ybar,ymax
         if(abs(ymax).gt.abs(ymin)) then
            pinl=kmax
         else
            pinl=kmin
         endif
         pinv=bufigram(pinl)
         dclevel=dble(ybar)
         frzpda=0.0d0
         sivcflag=0
         fvsi_calc=0.0d0 
c This is set to zero because there's no smoothed IGRAM for the AC case. 
c However, if you wish to calculate fvsi_calc value for the AC case, you
c should do it in the deglitching process.
         call vsub(bufigram,1,ybar,0,bufigram,1,nip)
         return
      else
         write(*,*)'DC Interferogram: ymin,ybar,ymax ',ymin,ybar,ymax
      endif

c  Smooth DC igram (unshifted with respect to the original)
      call smooth_igram(nip,margin,bufigram,bufsmoo,frsp,dclevel,frzpda,
     & pinv,pinl)

c       write(*,*) 'pinl,nip,frzpda=',pinl,nip,frzpda

c  Find the max/min values of the smoothed igram (BUFSMOO).
      smoo_min=bufsmoo(1)
      smoo_max=bufsmoo(1)
      do jip=2,nip
         if(bufsmoo(jip).lt.smoo_min) smoo_min=bufsmoo(jip)
         if(bufsmoo(jip).gt.smoo_max) smoo_max=bufsmoo(jip)
      enddo

c  Calculate average value and standard deviation of BUFSMOO to
c  approximate the SIA and SIS values, in case they do not exist
c  from the solar tracker

      smoo_std =0
c     write(42,*)bufsmoo
c     smoo_mean=sum(dble(bufsmoo))/nip
c     do jip=1,nip
c        smoo_std =smoo_std+(bufsmoo(jip)-smoo_mean)**2
c     enddo
c     smoo_std=sqrt(smoo_std/nip)
c     fvsi_calc=abs(smoo_std/smoo_mean)
c     write(*,*)'fvsi_calc,smoo_mean=',fvsi_calc,smoo_mean

      tot = 0.0d0
      do jip=1,nip
         tot = tot + bufsmoo(jip)
      enddo
      smoo_mean=sngl( tot/nip )
      tot = 0.0d0
      do jip=1,nip
         tot = tot + (bufsmoo(jip)-smoo_mean)**2
      enddo
      smoo_std=sngl( dsqrt(tot/nip) )
      fvsi_calc=abs(smoo_std/smoo_mean)
c     write(*,*)'fvsi_calc,smoo_mean=',fvsi_calc,smoo_mean

c  If FRSP <=0, then return the AC igram.
      if(frsp.le.0.0) then  ! Subtract smoothed igram from original
         call vsub(bufigram,1,bufsmoo,1,bufigram,1,nip)
         sivcflag=1
         return
      endif

c The following if-statement does nothing, since the AC interferogram
c does not get smoothed. But it has been tested for the DC case
c and both methods return similar fvsi_calc values.
c     if( ybar.lt.0.2*(ymax-ymin) ) then      ! Treat as AC
c        fvsi_calc=smoo_std*1.6/(ymax-ymin)   ! Assume ~80% ME at ZPD
c     else ! Treat as DC
c        fvsi_calc=smoo_std/smoo_mean
c     endif

c     write(*,*)'smoo_mean,smoo_std,ymax,ymin,fvsi_calc,fvsi_calcAC',
c    & smoo_mean,smoo_std,ymax,ymin,fvsi_calc,smoo_std*1.6/(ymax-ymin)

c  If the smoothed igram crosses zero, or even comes close to zero,
c  we don't want to divide anything by it, because this would amplify
c  the noise and could even cause an Inf or Nan.
c  So we add a small offset to the smoothed igram that is just sufficient to
c  to keep the smoothed igram the same sign, either entirely +ve or -ve.
      if(smoo_min/ybar.lt.sthresh) then
         soffset=sthresh*ybar-smoo_min
      elseif(smoo_max/ybar.lt.sthresh) then
         soffset=sthresh*ybar-smoo_max
      else
         soffset=0.0
      endif

c  Perform the SIV correction on the interferogram, and convert to AC.
      do jip=1,nip    !  Perform SIV correction 
         bufigram(jip)=sngl( dclevel*
     &   ( (soffset+bufigram(jip))/(soffset+bufsmoo(jip)) -1. ) )
      end do
      sivcflag=2
      return
      end
