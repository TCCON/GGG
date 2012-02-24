      subroutine resample(vin,foff,m1,vout,nmp)
c  Stretches the input vector VIN by a factor 1+FOFF and
c  then resamples it by Lagrangian interpolation
c  such that    VOUT[m1+i] = VIN[(M1+i)*(1+FOFF)]
c  or           VOUT[i] = VIN[(M1+i)*(1+FOFF)-M1]
c  or           VOUT[i] = VIN[i*(1+FOFF)+M1*FOFF]
c
c  Inputs:
c      VIN(NMP)   R*4  Input vector
c      FOFF       R*4  Fractional offset (usually ~ 0.000001 )
c      M1         I*4  Absolute index of VIN(0)
c      NMP        I*4  Number of Points
c
c Outputs:
c      VOUT(NMP)  R*4  Output Vector
c
c The length of the resampling operator is 2*NTERM+1. so we lose
c NTERM points from each end of the input vector, even when the
c frequency stretch is zero.
c
c  To avoid array bound violations in VIN, we must ensure that
c      1+NTERM =<  i*(1+FOFF)+M1*FOFF =< NMP-NTERM
c
c  This means that:
c      i >= (1+NTERM-M1*FOFF)/(1+FOFF)
c      i =< (NMP-NTERM-M1*FOFF)/(1+FOFF)
c
c  Hence
c      i_min = NINT((1+NTERM-M1*FOFF)/(1+FOFF))
c      i_max = NINT((NMP-NTERM-M1*FOFF)/(1+FOFF))
c
      implicit none
      integer*4 m1,nmp,i,is,i_min,i_max,nterm
      real*4 vin(nmp),vout(nmp),foff,shift,fs,polyinterp
c
      nterm=4
c      i_min=max0(1,1+nterm-nint(foff*(m1+1+nterm+1)))
c      i_max=min0(nmp,nmp-nterm-nint(foff*(m1+nmp-nterm-1)))
      i_min=max0(1,nint((1+nterm-m1*foff)/(1.+foff)))
      i_max=min0(nmp,nint((nmp-nterm-m1*foff)/(1.+foff)))
c      write(*,*)'resample',m1,i_min,i_max,foff,nmp,nint(foff*(m1+1))
c
c Interpolate to the shifted/stretched frequencies.
      do i=i_min,i_max
         shift=foff*(i+m1)
         is=nint(shift)
         fs=shift-is
         is=is+i-nterm
c         if(is.eq.0) write(*,*)'Resample:',m1,foff,shift,i_min,i,i_max
         vout(i)=polyinterp(vin(is),nterm,fs)
      end do

c  Initialize the first few points of the array which otherwise would all
c  be zero because of the shift or the width of the interpolation operator
      do i=i_min-1,1,-1
c         vout(i)=2*vout(i+1)-vout(i+2)
         vout(i)=vin(i)
      end do

c  Initialize the last few points of the array which otherwise would all
c  be zero because of the shift or the width of the interpolation operator
      do i=i_max+1,nmp
c         vout(i)=2*vout(i-1)-vout(i-2)
         vout(i)=vin(i)
      end do
c
      return
      end
