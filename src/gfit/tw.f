      real*4 vvv
      real*8 yyy, wlimit

      vvv=1600.0
      do k=1,12
         vvv=-1.5*vvv
         yyy=wlimit(vvv,'f9.4')
         write(*,*)k,vvv,yyy
      end do
      stop
      end

      include '../comn/lnbc.f'
      include 'wlimit.f'


c      function wlimit(vvv,wformat)
cc  Truncates the value of vvv so that it may be written
cc  in wformat without producing **** due to format overflow
cc
cc  Inputs:
cc     vvv      R*4  Floating point value
cc     wformat  C*(*)  Format (e.g. 'f5.3')
cc
cc  Output:
cc     wlimit   R*4  Truncated value of vvv
cc
cc  Examples:
cc    vvv=121.5,  wformat=f4.1', wlimit=99.9
cc    vvv=121.5,  wformat=f5.2', wlimit=121.5
cc    vvv=-21.5,  wformat=f5.2', wlimit=-9.99
cc    vvv=-21.5,  wformat=f5.1', wlimit=-21.5
cc
c
c      implicit none
c      integer*4 i,j,k1,k2,lw,lnbc
c      real*4 vvv
c      real*8 wmin,wmax,wlimit
c      character wformat*(*)
c
c      lw=lnbc(wformat)
cc
cc  Extract numerical values of k1 (=5) and k2 (=3) from wformat (='f5.3')
cc  This is *not* done by an internal read because the function wlimit is
cc  often used in write statements which would cause the G77 compiler to
cc  complain:  "I/O recursion: I/O started while already doing I/O"
cc  So reading K1 and K2 has to be done the old-fashioned way.
cc
cc  Skip the first character ('f') and then evaluate K1. 
c      k1=0
c      do i=2,lw
c        if(wformat(i:i).eq.'.') exit
c        k1=10*k1+ichar(wformat(i:i))-48
c      end do
c
cc  Evaluate K2
c      k2=0
c      do j=i+1,lw
c        k2=10*k2+ichar(wformat(j:j))-48
c      end do
c
cc Determine largest & smallest writable values allowed by wformat
c      wmax=dfloat(10**(k1-1)-1)/10**k2
c      wmin=-dfloat(10**(k1-2)-1)/10**k2
c
cc Truncate data value to wmax/wmin if necessary.
c      wlimit=vvv
c      if(wlimit.gt.wmax) wlimit=wmax
c      if(wlimit.lt.wmin) wlimit=wmin
c      write(*,*)wformat,vvv,k1,k2,wmin,wmax,wlimit
c      return
c      end
