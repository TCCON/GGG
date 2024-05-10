      subroutine resample_ifg(nip, bufigram,
     & bufigram_out,
     & evenodd, deltaX, halfSincSize)

c  Resamples the DC-corrected interferogram, given pre-determined 
c  shift (deltaX) and parity of points to shift (evenodd)
c 
c  Written for i2s by Vanessa Sherlock based on KIT algorithms from 
c  Frank Hase and Michael Gisi. Slightly modified by GCT.
c
c  Inputs:
c           nip           I*4  Number of interferogram points
c           bufigram(nip) R*4  Raw Interferogram (AC or DC)
c           evenodd       I*4  which sample points to shift 
c           deltaX        R*4  interferogram shift
c           halfSincSize  I*4  should be > 50
c
c Outputs: 
c           bufigram_out(nip) R*4  Resampled Interferogram
c
c
c
      implicit none
      integer*4
     & i, j,
     & halfSincSize,nip,evenodd

      real*4
     & deltaX,
     & bufigram(nip), ! Data storage and FFT bufigram
     & bufigram_out(nip) ! Data storage and FFT bufigram

      real*8 
     & sinc(2*halfSincSize + 1), Pi, agmt, wert, u

      Pi=DACOS(-1.D0) 

c calc. modified sinc
      do i=0,2*halfSincSize
         u=1.D0*(i - halfSincSize)
         agmt=ABS(Pi*(u + dble(deltaX))) + 1.e-10
c       write(*,*)'u,agmt=',u,agmt
         sinc(i+1) = sin(agmt)/agmt*0.5*(1. + 
     &                 cos(Pi*u/(1.d0*halfSincSize)))
c       write(*,*) i, sinc(i+1),agmt,sin(agmt)/agmt,deltaX
      end do

c resample
      do i=1,halfSincSize
         bufigram_out(i)=0.0
      enddo

      do i=nip-halfSincSize+1, nip
         bufigram_out(i)=0.0
      enddo

      do i=halfSincSize+1,nip-halfSincSize 
         if (MOD(i,2) .eq. evenodd) then
            bufigram_out(i)=bufigram(i)
         endif 
         if (MOD(i,2) .ne. evenodd) then
            wert=0.0d0
            do j=0, 2*halfSincSize
               wert=wert+bufigram(i+j-halfSincSize)*sinc(j+1)
            enddo
            bufigram_out(i)=sngl(wert) 
         endif 
      enddo

cc write resampled interferogram into i/o bufigram
c      do i=1,nip
c         bufigram(i)=bufigram_out(i)
c      enddo

      return
      end
