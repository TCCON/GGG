c modfft.f
c
c  Revision 1.7  2016/12/09
c  Replaced all instances of "24" with MFFT
c  and added MFFT to the common statement.
c  Now MFFT is defined in a single place.
c
c  Revision 1.6  2012/03/11
c  Replaced the obsolete arithmetic if-statements, e.g.
c      if (m-n8pow*3-1) 50, 40, 30                                       
c  to prevent compiler warnings such as:
c  Warning: Obsolete: arithmetic IF statement at (1)
c  This was done in 9 different places.
c  Results appear identical to those previously.
c
c  Also added STOP statments when NFFT exceeds 2**24.
c  Previously the program returned after writing a warning message

c
c  Revision 1.5  2004/05/14  JFB
c    Replaced leading tabs with spaces
c    Fixed indentation of variable declarations
c    Added "implicit none" to individual subroutines
c    Replaced two instances of "pii = 4.*datan(1.d0) with
c      "pii = real(4.d0*datan(1.d0))" to avoid ftnchek warnings
c    Removed unused label "90" in front of "call ord1 (m, b)"
c    Added parenthesis around "-p7", "-tr3", and "-t8" to avoid G77 warnings
c
c  Revision 1.4  2002/09/16  GCT
c    Extended code to support 2^24 transforms.
c
c  Revision 1.3  1998/11/05  GCT
c    Replaced most of the arrays dimensioned (*) by the actual adjustable
c    parameters so that the array bound checker can report any violations.
c    Also replace * by & as the continuation character.
c    Also defined pii = 4*datan(1.0d0) instead of 4*atan(1.)
c    On the Sun the latter was only good to 6 sig figs (3.14159).
c
c  Revision 1.2  1996/10/20  GCT
c   -Extended code to support 2^22 transforms for MkIV 120 cm InSb runs
c   -Disabled 'feature' which returned DC term in B(NFFT+1) and the
c    Nyquist term in B(NFFT+2), requiring declared dimension of B(NFFT+2).
c    Routine now requires declared dimension of B(NFFT).
c   -Removed do loops, introduced by MCA, to perform origin shift by
c    multiplying every other complex spectral point by -1.
c   -Routine I/O now corresponds exactly to original routine with
c    the array packed in the conventional (wrapped around) FFT format.
c
c Revision 1.1  1996/04/30  23:36:20  aychang
c S0146  clean version.
c        added implicit delcarations to compile with type checking
c
c Revision 1.0  1996/01/26  00:42:47  aychang
c Initial revision
c
c  L. Delbouille, July 26, 1987; modified by JWB 4/13/91; MCA 5/11/92
c
c  The original author of the FFT subroutines below is: G. D. Bergland,
c  "A Radix-Eight Fast Fourier Transform Subroutine for Real-Valued Series",
c  IEEE Trans. Audio Electroacoust., vol. AU-17, pp 138-144, June 1969
c
c ------------------------------------------------------------------
c Subroutine:  ffak              Fast fourier analysis subroutine
c-------------------------------------------------------------------
      subroutine ffak (b, nfft)

c  This subroutine replaces the real vector b(k), (k=1,2,...,n), with 
c  its finite discrete fourier transform.  The dc term is returned in 
c  location b(1).  Thereafter, the j'th harmonic is
c  returned as a complex number stored as b(2*j+1) + i b(2*j+2).
c  Note that the n/2 harmonic (Nyquist) is returned in b(2). 
c  The subroutine is called as ffak (b,n) where n=2**m and b is an
c  n term real array.  A real-valued, radix 8 algorithm is used with
c  in-place reordering and the trig functions are computed as needed.

      implicit none
      integer*4 nfft
      real*4 b(nfft)

      integer n, i, m, nn, int, n8pow, it, mfftpow
      parameter (mfftpow=24)

      real*4 pi8

      real*4 cos, sin

      real*4 pii, p7, p7two, c22, s22, pi2
      common /con/ pii, p7, p7two, c22, s22, pi2

      pii = 4*atan(1.)
      pi8 = pii/8.
      p7 = 1./sqrt(2.)
      p7two = 2.*p7
      c22 = cos(pi8)
      s22 = sin(pi8)
      pi2 = 2.*pii
      n = 1
      do i=1,mfftpow
         m = i
         n = n*2
         if (n.eq.nfft) go to 20
      enddo
      stop 'ffak:  nfft not a power of 2 <= 2^24 '
c      return

  20  n8pow = m/3

c  do a radix 2 or radix 4 iteration first if one is required
      if (m.gt.n8pow*3+1) then        !  GCT 201203122
         nn = 4
         int = n/nn
         call r4tr (int, b(1), b(int+1), b(2*int+1), b(3*int+1))
      elseif(m.eq.n8pow*3+1) then
         nn = 2
         int = n/nn
         call r2tr (int, b(1), b(int+1))
      else
         nn = 1
      endif

c   Perform radix 8 iterations
      if (n8pow .gt. 0) then
         do it=1,n8pow
            nn = nn*8
            int = n/nn
            call r8trk (int,nn,b(1),b(int+1),b(2*int+1),b(3*int+1),
     &        b(4*int+1), b(5*int+1), b(6*int+1), b(7*int+1), b(1),
     &        b(int+1),b(2*int+1),b(3*int+1),b(4*int+1),b(5*int+1),
     &        b(6*int+1), b(7*int+1))
         enddo
      end if
      ! perform in-place reordering
      call ord1 (m, b)
      call ord2k (m, b)
      return
      end
c
c-----------------------------------------------------------------------
c subroutine:  r2tr      ^I radix 2 iteration subroutine
c-----------------------------------------------------------------------
      subroutine r2tr(int, b0, b1)

      implicit none
      integer*4 int
      real*4 b0(int), b1(int)

      integer*4 k
      real*4 t

      do k=1,int
         t = b0(k) + b1(k)
         b1(k) = b0(k) - b1(k)
         b0(k) = t
      enddo
      return
      end
c-----------------------------------------------------------------------
c subroutine:  r4tr                 radix 4 iteration subroutine
c-----------------------------------------------------------------------
      subroutine r4tr (int, b0, b1, b2, b3)

      implicit none
      integer*4 int
      real*4 b0(int), b1(int), b2(int), b3(int)

      integer*4 k
      real*4 r0, r1

      do k=1,int
         r0 = b0(k) + b2(k)
         r1 = b1(k) + b3(k)
         b2(k) = b0(k) - b2(k)
         b3(k) = b1(k) - b3(k)
         b0(k) = r0 + r1
         b1(k) = r0 - r1
      end do
      return
      end
c-----------------------------------------------------------------------
c subroutine: r8trk                   radix 8 iteration subroutine
c for a number of points greater than 32k, up to 2048k
c-----------------------------------------------------------------------
      subroutine r8trk (int, nn, br0, br1, br2, br3, br4, br5, br6,
     &    br7, bi0, bi1, bi2, bi3, bi4, bi5, bi6, bi7)

      implicit none
      integer*4 int, nn, mfftpow
      parameter (mfftpow=24)
      real*4 br0(*), br1(*), br2(*), br3(*), br4(*),
     &  br5(*), br6(*), br7(*), bi0(*), bi1(*), bi2(*),
     &  bi3(*), bi4(*), bi5(*), bi6(*), bi7(*)

      integer*4 l(mfftpow), l1, l2, l3, l4, l5, l6, l7, l8, l9, l10,
     &  l11, l12, l13, l14, l15, l16, l17, l18, l19, l20,
     &  l21, l22, l23, l24

      integer*4 j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11, j12, 
     &  j13, j14, j15, j16, j17, j18, j19, j20, j21, j22, j23,
     &  jthet, k, ji, jr, jl, k0, kl, int8, jlast, j, j0, th2

      real*4 t0, t1, t2, t3, t4, t5, t6, t7, tr0, ti0,
     &  tr1, ti1, tr2, ti2, tr3, ti3, tr4, ti4, tr5, ti5,
     &  tr6, ti6, tr7, ti7, c1, c2, c3, c4, c5, c6, c7,
     &  s1, s2, s3, s4, s5, s6, s7, piovn, pr, pi, arg

      real*4 sin, cos

      real*4 pii, p7, p7two, c22, s22, pi2
      common /con/ pii, p7, p7two, c22, s22, pi2

      equivalence
     & (l24,l(1)),
     & (l23,l(2)),
     & (l22,l(3)),
     & (l21,l(4)),
     & (l20,l(5)),
     & (l19,l(6)),
     & (l18,l(7)),
     & (l17,l(8)),
     & (l16,l(9)),
     & (l15,l(10)),
     & (l14,l(11)),
     & (l13,l(12)),
     & (l12,l(13)),
     & (l11,l(14)),
     & (l10,l(15)),
     & (l9,l(16)),
     & (l8,l(17)),
     & (l7,l(18)),
     & (l6,l(19)),
     & (l5,l(20)),
     & (l4,l(21)),
     & (l3,l(22)),
     & (l2,l(23)),
     & (l1,l(24))

c  set up counters such that jthet steps through the arguments of w,
c  jr steps through starting locations for the real part of the
c  intermediate results and ji steps through starting locations
c  of the imaginary part of the intermediate results.

      l(1) = nn/8
      do k=2,mfftpow
         if (l(k-1).gt.2) then        !  GCT 201203122
            l(k) = l(k-1)/2
         else
            if (l(k-1).lt.2) l(k-1) = 2
            l(k) = 2
         endif
      end do
      piovn = pii/float(nn)
      ji = 3
      jl = 2
      jr = 2
      do 120 j1=2,l1,2
      do 120 j2=j1,l2,l1
      do 120 j3=j2,l3,l2
      do 120 j4=j3,l4,l3
      do 120 j5=j4,l5,l4
      do 120 j6=j5,l6,l5
      do 120 j7=j6,l7,l6
      do 120 j8=j7,l8,l7
      do 120 j9=j8,l9,l8
      do 120 j10=j9,l10,l9
      do 120 j11=j10,l11,l10
      do 120 j12=j11,l12,l11
      do 120 j13=j12,l13,l12
      do 120 j14=j13,l14,l13
      do 120 j15=j14,l15,l14
      do 120 j16=j15,l16,l15
      do 120 j17=j16,l17,l16
      do 120 j18=j17,l18,l17
      do 120 j19=j18,l19,l18
      do 120 j20=j19,l20,l19
      do 120 j21=j20,l21,l20
      do 120 j22=j21,l22,l21
      do 120 j23=j22,l23,l22
      do 120 jthet=j23,l24,l23
         th2 = jthet - 2
         if (th2 .gt. 0) go to 90
         do k=1,int
            t0 = br0(k) + br4(k)
            t1 = br1(k) + br5(k)
            t2 = br2(k) + br6(k)
            t3 = br3(k) + br7(k)
            t4 = br0(k) - br4(k)
            t5 = br1(k) - br5(k)
            t6 = br2(k) - br6(k)
            t7 = br3(k) - br7(k)
            br2(k) = t0 - t2
            br3(k) = t1 - t3
            t0 = t0 + t2
            t1 = t1 + t3
            br0(k) = t0 + t1
            br1(k) = t0 - t1
            pr = p7*(t5-t7)
            pi = p7*(t5+t7)
            br4(k) = t4 + pr
            br7(k) = t6 + pi
            br6(k) = t4 - pr
            br5(k) = pi - t6
         end do
         if (nn-8 .le. 0) go to 120

         k0 = int*8 + 1
         kl = k0 + int - 1
         do k=k0,kl
            pr = p7*(bi2(k)-bi6(k))
            pi = p7*(bi2(k)+bi6(k))
            tr0 = bi0(k) + pr
            ti0 = bi4(k) + pi
            tr2 = bi0(k) - pr
            ti2 = bi4(k) - pi
            pr = p7*(bi3(k)-bi7(k))
            pi = p7*(bi3(k)+bi7(k))
            tr1 = bi1(k) + pr
            ti1 = bi5(k) + pi
            tr3 = bi1(k) - pr
            ti3 = bi5(k) - pi
            pr = tr1*c22 - ti1*s22
            pi = ti1*c22 + tr1*s22
            bi0(k) = tr0 + pr
            bi6(k) = tr0 - pr
            bi7(k) = ti0 + pi
            bi1(k) = pi - ti0
            pr = (-tr3)*s22 - ti3*c22
            pi = tr3*c22 - ti3*s22
            bi2(k) = tr2 + pr
            bi4(k) = tr2 - pr
            bi5(k) = ti2 + pi
            bi3(k) = pi - ti2
         end do
         go to 120

  90     arg = float(th2)*piovn
         c1 = cos(arg)
         s1 = sin(arg)
         c2 = c1**2 - s1**2
         s2 = c1*s1 + c1*s1
         c3 = c1*c2 - s1*s2
         s3 = c2*s1 + s2*c1
         c4 = c2**2 - s2**2
         s4 = c2*s2 + c2*s2
         c5 = c2*c3 - s2*s3
         s5 = c3*s2 + s3*c2
         c6 = c3**2 - s3**2
         s6 = c3*s3 + c3*s3
         c7 = c3*c4 - s3*s4
         s7 = c4*s3 + s4*c3
         int8 = int*8
         j0 = jr*int8 + 1
         k0 = ji*int8 + 1
         jlast = j0 + int - 1
         do j=j0,jlast
            k = k0 + j - j0
            tr1 = br1(j)*c1 - bi1(k)*s1
            ti1 = br1(j)*s1 + bi1(k)*c1
            tr2 = br2(j)*c2 - bi2(k)*s2
            ti2 = br2(j)*s2 + bi2(k)*c2
            tr3 = br3(j)*c3 - bi3(k)*s3
            ti3 = br3(j)*s3 + bi3(k)*c3
            tr4 = br4(j)*c4 - bi4(k)*s4
            ti4 = br4(j)*s4 + bi4(k)*c4
            tr5 = br5(j)*c5 - bi5(k)*s5
            ti5 = br5(j)*s5 + bi5(k)*c5
            tr6 = br6(j)*c6 - bi6(k)*s6
            ti6 = br6(j)*s6 + bi6(k)*c6
            tr7 = br7(j)*c7 - bi7(k)*s7
            ti7 = br7(j)*s7 + bi7(k)*c7

            t0 = br0(j) + tr4
            t1 = bi0(k) + ti4
            tr4 = br0(j) - tr4
            ti4 = bi0(k) - ti4
            t2 = tr1 + tr5
            t3 = ti1 + ti5
            tr5 = tr1 - tr5
            ti5 = ti1 - ti5
            t4 = tr2 + tr6
            t5 = ti2 + ti6
            tr6 = tr2 - tr6
            ti6 = ti2 - ti6
            t6 = tr3 + tr7
            t7 = ti3 + ti7
            tr7 = tr3 - tr7
            ti7 = ti3 - ti7

            tr0 = t0 + t4
            ti0 = t1 + t5
            tr2 = t0 - t4
            ti2 = t1 - t5
            tr1 = t2 + t6
            ti1 = t3 + t7
            tr3 = t2 - t6
            ti3 = t3 - t7
            t0 = tr4 - ti6
            t1 = ti4 + tr6
            t4 = tr4 + ti6
            t5 = ti4 - tr6
            t2 = tr5 - ti7
            t3 = ti5 + tr7
            t6 = tr5 + ti7
            t7 = ti5 - tr7
            br0(j) = tr0 + tr1
            bi7(k) = ti0 + ti1
            bi6(k) = tr0 - tr1
            br1(j) = ti1 - ti0
            br2(j) = tr2 - ti3
            bi5(k) = ti2 + tr3
            bi4(k) = tr2 + ti3
            br3(j) = tr3 - ti2
            pr = p7*(t2-t3)
            pi = p7*(t2+t3)
            br4(j) = t0 + pr
            bi3(k) = t1 + pi
            bi2(k) = t0 - pr
            br5(j) = pi - t1
            pr = (-p7)*(t6+t7)
            pi = p7*(t6-t7)
            br6(j) = t4 + pr
            bi1(k) = t5 + pi
            bi0(k) = t4 - pr
            br7(j) = pi - t5
         end do
         jr = jr + 2
         ji = ji - 2
         if (ji-jl .gt. 0) go to 120
         ji = 2*jr - 1
         jl = jr
 120  continue

      return
      end
c-----------------------------------------------------------------------
c subroutine:  ord1
c-----------------------------------------------------------------------
      subroutine ord1 (m,b)

      implicit none
      integer*4 m
      real*4 b(2**m)

      integer*4 k, kl, j, n
      real*4 t

      k = 4
      kl = 2
      n = 2**m
      do j=4,n,2
         if (k.gt.j) then        !  GCT 20120311
            t = b(j)
            b(j) = b(k)
            b(k) = t
         endif       ! GCT 20120311
         k = k - 2
         if (k.gt.kl) cycle        !  GCT 20120311
         k = 2*j
         kl = j
      end do
      return
      end
c-----------------------------------------------------------------------
c subroutine:  ord2k
c in-place reordering subroutine, for up to 2 million points
c-----------------------------------------------------------------------
      subroutine ord2k (m, b)

      implicit none
      integer*4 m,mfftpow
      parameter (mfftpow=24)
      real*4 b(2**m)

      integer*4 l(mfftpow), l1, l2, l3, l4, l5, l6, l7, l8, l9, l10,
     &  l11, l12, l13, l14, l15, l16, l17, l18, l19, l20,
     &  l21, l22, l23, l24

      integer*4 j1, j2, j3, j4, j5, j6, j7, j8, j9, j10,
     &  j11, j12, j13, j14, j15, j16, j17, j18, j19, j20,
     &  j21, j22, j23, n, k, ji, ij

      real*4 t

      equivalence
     & (l24,l(1)),
     & (l23,l(2)),
     & (l22,l(3)),
     & (l21,l(4)),
     & (l20,l(5)),
     & (l19,l(6)),
     & (l18,l(7)),
     & (l17,l(8)),
     & (l16,l(9)),
     & (l15,l(10)),
     & (l14,l(11)),
     & (l13,l(12)),
     & (l12,l(13)),
     & (l11,l(14)),
     & (l10,l(15)),
     & (l9,l(16)),
     & (l8,l(17)),
     & (l7,l(18)),
     & (l6,l(19)),
     & (l5,l(20)),
     & (l4,l(21)),
     & (l3,l(22)),
     & (l2,l(23)),
     & (l1,l(24))

      n = 2**m
      l(1) = n
      do k=2,m
         l(k) = l(k-1)/2
      enddo
      do k=m,23
         l(k+1) = 2
      enddo
      ij = 2
      do 40 j1=2,l1,2
      do 40 j2=j1,l2,l1
      do 40 j3=j2,l3,l2
      do 40 j4=j3,l4,l3
      do 40 j5=j4,l5,l4
      do 40 j6=j5,l6,l5
      do 40 j7=j6,l7,l6
      do 40 j8=j7,l8,l7
      do 40 j9=j8,l9,l8
      do 40 j10=j9,l10,l9
      do 40 j11=j10,l11,l10
      do 40 j12=j11,l12,l11
      do 40 j13=j12,l13,l12
      do 40 j14=j13,l14,l13
      do 40 j15=j14,l15,l14
      do 40 j16=j15,l16,l15
      do 40 j17=j16,l17,l16
      do 40 j18=j17,l18,l17
      do 40 j19=j18,l19,l18
      do 40 j20=j19,l20,l19
      do 40 j21=j20,l21,l20
      do 40 j22=j21,l22,l21
      do 40 j23=j22,l23,l22
      do 40 ji=j23,l24,l23
         if (ij.lt.ji) then        !  GCT 20120311
            t = b(ij-1)
            b(ij-1) = b(ji-1)
            b(ji-1) = t
            t = b(ij)
            b(ij) = b(ji)
            b(ji) = t
         endif                     ! GCT 20120311
  40     ij = ij + 2
      return
      end

c-----------------------------------------------------------------------
c subroutine:  ffsk
c fast fourier synthesis subroutine, radix 8-4-2
c  modified by JWB 8/10/91 for 2**21 points
c-----------------------------------------------------------------------
      subroutine ffsk(b, nfft)

c This subroutine synthesizes the real vector b(k), where k=1,2,...,n. 
c The initial fourier coefficients are placed in the b array of size n.  
c The jth harmonic is stored as b(2*j+1) + i b(2*j+2).
c The dc term is in b(1)
c The n/2 harmonic is in b(2).
c The subroutine is called as ffsk(b,n) where n=2**m and b is the n term 
c real array discussed above.

      implicit none
      integer*4 nfft,mfftpow
      parameter (mfftpow=24)
      real*4 b(nfft)

      integer*4 i, n, nn, n8pow, it, int, m
      real*4 con, pi8

      real*4 cos, sin

      real*4 pii, p7, p7two, c22, s22, pi2
      common /con1/ pii, p7, p7two, c22, s22, pi2

      pii = 4*atan(1.)
      pi8 = pii/8.
      p7 = 1./sqrt(2.)
      p7two = 2.*p7
      c22 = cos(pi8)
      s22 = sin(pi8)
      pi2 = 2.*pii
      n = 1
      do i=1,mfftpow
         m = i
         n = n*2
         if (n.eq.nfft) go to 20
      enddo
      write(*,*)'mfftpow, 2**mfftpow, nfft=',mfftpow,2**mfftpow,nfft
      stop 'ffsk: nfft not a power of 2 <= 2^24 '

20    continue
      con = 1.0/float(nfft)
      do i=1,nfft
         b(i) = b(i)*con
      enddo

      n8pow = m/3

c reorder the input Fourier coefficients
      call ord2k(m, b)
      call ord1(m, b)


c perform the radix 8 iterations
      if (n8pow .gt. 0) then
         nn = n
         do it=1,n8pow
c            write (*,'(a1,i1,$)') '^H',it
            int = n/nn
            call r8synk(int, nn, b, b(int+1), b(2*int+1), b(3*int+1),
     &        b(4*int+1), b(5*int+1), b(6*int+1), b(7*int+1), b(1),
     &        b(int+1), b(2*int+1), b(3*int+1), b(4*int+1), b(5*int+1),
     &        b(6*int+1), b(7*int+1))
            nn = nn/8
         enddo

      end if

c  Do a radix 2 or radix 4 iteration if one is required
      if (m.gt.n8pow*3+1) then        ! GCT 20120311
         int = n/4
         call r4syn(int, b(1), b(int+1), b(2*int+1), b(3*int+1))
      elseif (m.eq.n8pow*3+1) then
         int = n/2
         call r2tr(int, b(1), b(int+1))
      endif
      return
      end
c
c-----------------------------------------------------------------------
c subroutine:  r8synk
c radix 8 synthesis subroutine
c-----------------------------------------------------------------------
      subroutine r8synk(int, nn, br0,br1,br2,br3,br4,br5,br6,br7,
     &    bi0, bi1, bi2, bi3, bi4, bi5, bi6, bi7)

      implicit none
      integer*4 int, nn, mfftpow
      parameter(mfftpow=24)
      real*4 br0(*), br1(*), br2(*), br3(*), br4(*),
     &  br5(*), br6(*), br7(*), bi0(*), bi1(*), bi2(*),
     &  bi3(*), bi4(*), bi5(*), bi6(*), bi7(*)

      integer*4 l(mfftpow), l1, l2, l3, l4, l5, l6, l7, l8, l9, l10,
     &  l11, l12, l13, l14, l15, l16, l17, l18, l19, l20,
     &  l21, l22, l23, l24

      integer*4 j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11, j12,
     &  j13, j14, j15, j16, j17, j18, j19, j20, j21, j22, j23,
     &  jthet, k, ji, jr, jl, k0, kl, int8, jlast, j, j0, th2 

      real*4 t0, t1, t2, t3, t4, t5, t6, t7, t8, tr0, ti0,
     &  tr1, ti1, tr2, ti2, tr3, ti3, tr4, ti4, tr5, ti5,
     &  tr6, ti6, tr7, ti7, c1, c2, c3, c4, c5, c6, c7,
     &  s1, s2, s3, s4, s5, s6, s7,
     &  piovn, pr, pi, tt0, tt1, rr, ri, ttr6, ttr7, arg

      real*4 sin, cos

      real*4 pii, p7, p7two, c22, s22, pi2
      common /con1/ pii, p7, p7two, c22, s22, pi2

      equivalence
     & (l24,l(1)),
     & (l23,l(2)),
     & (l22,l(3)),
     & (l21,l(4)),
     & (l20,l(5)),
     & (l19,l(6)),
     & (l18,l(7)),
     & (l17,l(8)),
     & (l16,l(9)),
     & (l15,l(10)),
     & (l14,l(11)),
     & (l13,l(12)),
     & (l12,l(13)),
     & (l11,l(14)),
     & (l10,l(15)),
     & (l9,l(16)),
     & (l8,l(17)),
     & (l7,l(18)),
     & (l6,l(19)),
     & (l5,l(20)),
     & (l4,l(21)),
     & (l3,l(22)),
     & (l2,l(23)),
     & (l1,l(24))

      l(1) = nn/8
      do k=2,mfftpow
         if (l(k-1).le.2) then      !  GCT 20120311
            if(l(k-1).lt.2) l(k-1) = 2
            l(k) = 2
         else
            l(k) = l(k-1)/2
         endif
      end do
      piovn = pii/float(nn)
      ji = 3
      jl = 2
      jr = 2

      do 120 j1=2,l1,2
      do 120 j2=j1,l2,l1
      do 120 j3=j2,l3,l2
      do 120 j4=j3,l4,l3
      do 120 j5=j4,l5,l4
      do 120 j6=j5,l6,l5
      do 120 j7=j6,l7,l6
      do 120 j8=j7,l8,l7
      do 120 j9=j8,l9,l8
      do 120 j10=j9,l10,l9
      do 120 j11=j10,l11,l10
      do 120 j12=j11,l12,l11
      do 120 j13=j12,l13,l12
      do 120 j14=j13,l14,l13
      do 120 j15=j14,l15,l14
      do 120 j16=j15,l16,l15
      do 120 j17=j16,l17,l16
      do 120 j18=j17,l18,l17
      do 120 j19=j18,l19,l18
      do 120 j20=j19,l20,l19
      do 120 j21=j20,l21,l20
      do 120 j22=j21,l22,l21
      do 120 j23=j22,l23,l22
      do 120 jthet=j23,l24,l23
      th2 = jthet - 2
      if (th2.gt.0.0) go to 90        !  GCT 201203122
      do k=1,int
         t0 = br0(k) + br1(k)
         t1 = br0(k) - br1(k)
         t2 = br2(k) + br2(k)
         t3 = br3(k) + br3(k)
         t4 = br4(k) + br6(k)
         t6 = br7(k) - br5(k)
         t5 = br4(k) - br6(k)
         t7 = br7(k) + br5(k)
         pr = p7*(t7+t5)
         pi = p7*(t7-t5)
         tt0 = t0 + t2
         tt1 = t1 + t3
         t2 = t0 - t2
         t3 = t1 - t3
         t4 = t4 + t4
         t5 = pr + pr
         t6 = t6 + t6
         t7 = pi + pi
         br0(k) = tt0 + t4
         br1(k) = tt1 + t5
         br2(k) = t2 + t6
         br3(k) = t3 + t7
         br4(k) = tt0 - t4
         br5(k) = tt1 - t5
         br6(k) = t2 - t6
         br7(k) = t3 - t7
      end do
      if (nn.le.8) cycle        !  GCT 201203122
      k0 = int*8 + 1
      kl = k0 + int - 1
      do k=k0,kl
         t1 = bi0(k) + bi6(k)
         t2 = bi7(k) - bi1(k)
         t3 = bi0(k) - bi6(k)
         t4 = bi7(k) + bi1(k)
         pr = t3*c22 + t4*s22
         pi = t4*c22 - t3*s22
         t5 = bi2(k) + bi4(k)
         t6 = bi5(k) - bi3(k)
         t7 = bi2(k) - bi4(k)
         t8 = bi5(k) + bi3(k)
         rr = t8*c22 - t7*s22
         ri = (-t8)*s22 - t7*c22
         bi0(k) = (t1+t5) + (t1+t5)
         bi4(k) = (t2+t6) + (t2+t6)
         bi1(k) = (pr+rr) + (pr+rr)
         bi5(k) = (pi+ri) + (pi+ri)
         t5 = t1 - t5
         t6 = t2 - t6
         bi2(k) = p7two*(t6+t5)
         bi6(k) = p7two*(t6-t5)
         rr = pr - rr
         ri = pi - ri
         bi3(k) = p7two*(ri+rr)
         bi7(k) = p7two*(ri-rr)
      end do
      go to 120
  90  arg = float(th2)*piovn
      c1 = cos(arg)
      s1 = -sin(arg)
      c2 = c1**2 - s1**2
      s2 = c1*s1 + c1*s1
      c3 = c1*c2 - s1*s2
      s3 = c2*s1 + s2*c1
      c4 = c2**2 - s2**2
      s4 = c2*s2 + c2*s2
      c5 = c2*c3 - s2*s3
      s5 = c3*s2 + s3*c2
      c6 = c3**2 - s3**2
      s6 = c3*s3 + c3*s3
      c7 = c3*c4 - s3*s4
      s7 = c4*s3 + s4*c3
      int8 = int*8
      j0 = jr*int8 + 1
      k0 = ji*int8 + 1
      jlast = j0 + int - 1
      do j=j0,jlast
         k = k0 + j - j0
         tr0 = br0(j) + bi6(k)
         ti0 = bi7(k) - br1(j)
         tr1 = br0(j) - bi6(k)
         ti1 = bi7(k) + br1(j)
         tr2 = br2(j) + bi4(k)
         ti2 = bi5(k) - br3(j)
         tr3 = bi5(k) + br3(j)
         ti3 = bi4(k) - br2(j)
         tr4 = br4(j) + bi2(k)
         ti4 = bi3(k) - br5(j)
         t0  = br4(j) - bi2(k)
         t1  = bi3(k) + br5(j)
         tr5 = p7*(t0+t1)
         ti5 = p7*(t1-t0)
         tr6 = br6(j) + bi0(k)
         ti6 = bi1(k) - br7(j)
         t0  = br6(j) - bi0(k)
         t1  = bi1(k) + br7(j)
         tr7 = (-p7)*(t0-t1)
         ti7 = (-p7)*(t1+t0)
         t0  = tr0 + tr2
         t1  = ti0 + ti2
         t2  = tr1 + tr3
         t3  = ti1 + ti3
         tr2 = tr0 - tr2
         ti2 = ti0 - ti2
         tr3 = tr1 - tr3
         ti3 = ti1 - ti3
         t4  = tr4 + tr6
         t5  = ti4 + ti6
         t6  = tr5 + tr7
         t7  = ti5 + ti7
         ttr6 = ti4 - ti6
         ti6 = tr6 - tr4
         ttr7 = ti5 - ti7
         ti7 = tr7 - tr5
         br0(j) = t0 + t4
         bi0(k) = t1 + t5
         br1(j) = c1*(t2+t6) - s1*(t3+t7)
         bi1(k) = c1*(t3+t7) + s1*(t2+t6)
         br2(j) = c2*(tr2+ttr6) - s2*(ti2+ti6)
         bi2(k) = c2*(ti2+ti6) + s2*(tr2+ttr6)
         br3(j) = c3*(tr3+ttr7) - s3*(ti3+ti7)
         bi3(k) = c3*(ti3+ti7) + s3*(tr3+ttr7)
         br4(j) = c4*(t0-t4) - s4*(t1-t5)
         bi4(k) = c4*(t1-t5) + s4*(t0-t4)
         br5(j) = c5*(t2-t6) - s5*(t3-t7)
         bi5(k) = c5*(t3-t7) + s5*(t2-t6)
         br6(j) = c6*(tr2-ttr6) - s6*(ti2-ti6)
         bi6(k) = c6*(ti2-ti6) + s6*(tr2-ttr6)
         br7(j) = c7*(tr3-ttr7) - s7*(ti3-ti7)
         bi7(k) = c7*(ti3-ti7) + s7*(tr3-ttr7)
      end do
      jr = jr + 2
      ji = ji - 2
      if (ji.gt.jl) cycle        !  GCT 201203122
      ji = 2*jr - 1
      jl = jr
 120  continue
      return
      end

c-----------------------------------------------------------------------
c subroutine:  r4syn
c radix 4 synthesis
c-----------------------------------------------------------------------
      subroutine r4syn(int, b0, b1, b2, b3)

      implicit none
      integer*4 int
      real*4 b0(int), b1(int), b2(int), b3(int)

      integer*4 k
      real*4 t0, t1, t2, t3

      do k=1,int
         t0 = b0(k) + b1(k)
         t1 = b0(k) - b1(k)
         t2 = b2(k) + b2(k)
         t3 = b3(k) + b3(k)
         b0(k) = t0 + t2
         b2(k) = t0 - t2
         b1(k) = t1 + t3
         b3(k) = t1 - t3
      end do
      return
      end
c-----------------------------------------------------------------------
