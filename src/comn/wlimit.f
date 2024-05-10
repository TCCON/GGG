      function wlimit(vvv,wformat)
c  Returns the closest numerical value to vvv that can
c  be written as wformat (e.g. 'f9.4') without overflow.
c
c  Used to guarantee that a value written in f-format does
c  not overflow, producing '********' in the output file.
c
c  Inputs:
c     vvv      R*8  Floating point value
c     wformat  C*(*)  Format (e.g. 'f5.3')
c
c  Output:
c     wlimit   R**  Truncated value of vvv
c
c  Examples:
c    vvv=121.5,  wformat=f4.1', wlimit=99.9
c    vvv=121.5,  wformat=f5.2', wlimit=121.5
c    vvv=-21.5,  wformat=f5.2', wlimit=-9.99
c    vvv=-21.5,  wformat=f5.1', wlimit=-21.5
c    vvv=12345.  wformat=f9.4', wlimit=9999.9999
c
c  wlimit has to be double precision to handle cases where
c  the  output contains more than 7 significant digits.
c  Otherwise  9999.9999 would be rounded to 10000.0000
c  which would exceed the f9.4 format.


      implicit none
      integer*4 i,j,k1,k2,lw,lnblnk
      real*8 vvv
      real*8 wmin,wmax,wlimit
      character wformat*(*)

      lw=lnblnk(wformat)
c
c  Extract numerical values of k1 (=5) and k2 (=3) from wformat (='f5.3')
c  This is *not* done by an internal read because the function wlimit is
c  often used in write statements which would cause the G77 compiler to
c  complain:  "I/O recursion: I/O started while already doing I/O"
c  So reading K1 and K2 has to be done the old-fashioned way.
c
c  Skip the first character ('f') and then evaluate K1. 
      k1=0
      do i=2,lw
         if(wformat(i:i).eq.'.') exit
         k1=10*k1+ichar(wformat(i:i))-48
      end do

c  Evaluate K2
      k2=0
      do j=i+1,lw
         k2=10*k2+ichar(wformat(j:j))-48
      end do

c Determine largest & smallest writable values allowed by wformat
      wmax=dfloat(10**(k1-1)-1)/10**k2
      wmin=-dfloat(10**(k1-2)-1)/10**k2

c Truncate data value to wmax/wmin if necessary.
      wlimit=vvv
      if(wlimit.gt.wmax) wlimit=wmax
      if(wlimit.lt.wmin) wlimit=wmin
c      write(*,*)wformat,vvv,k1,k2,wmin,wmax,wlimit
      return
      end
