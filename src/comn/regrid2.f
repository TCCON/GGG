      subroutine regrid2(jinf,nin,fin,nhw,oper,odec,rdec,
     & koutf,nout,fout,nexpl,nexpr)
c
c  Convolves a uniformly-spaced input vector FIN with operator OPER.
c  Output FOUT is sampled onto a new uniform grid of spacing RDEC
c  times that of FIN. 
c
c  If OPER is a (truncated & windowed) unit-width SINC function,
c  then the result is simply FIN re-sampled onto the new RDEC grid.
c
c  INPUTS:
c     JINF       I*4  Absolute index of first point of FIN
c     NIN        I*4  Number of elements of FIN that are available
c     FIN(NIN)   R*4  Input vector (e.g. infinite-resolution spectrum)
c     NHW        I*4  Half-width of operator OPER wrt FIN
c     OPER(*)    R*4  Operator to be convolved with FIN
c     ODEC       I*4  Integer over-sampling factor of OPER wrt FIN
c     RDEC       R*8  Spacing of points in FOUT wrt FIN.
c     KOUTF      I*4  Absolute index of first point of FOUT
c     NOUT       I*4  Number of points in output vector
c
c  OUTPUTS:
c     FOUT(NOUT) R*4  Output Vector (Result of convolution)
c     NEXPL      I*4  Number of Extrapolated Points on Left 
c     NEXPR      I*4  Number of Extrapolated Points on Right
c------------------------------------------------------------------------
c
c    FOUT(koutf+i) = SUM_j {OPER[(koutf+i).RDEC-j].FIN(j)} i=1,Nout
c
c  Example 1: FIN is an infinite-resolution calculated transmittance
c  spectrum on a sub-doppler grid. OPER is the ILS of the measured
c  spectrum. FOUT will then be the calculated transmittance on the
c  grid of the measured spectrum.
c
c  Example 2: OPER contains a (truncated & windowed) SINC function,
c  narrower than any features in FIN. FOUT will be FIN resampled
c  onto a grid of spacing RDEC.
c  
c  This routine performs SINC resampling at the two sub-grid
c  points straddling each desired x-ordinate: xx=rdec*(k+kout)
c     x1 = int(xx*odec)/odec 
c     x2 = int(xx*odec+1)/odec 
c  And then linearly interpolates to x.
c     Y  = FIN(x1).(x2-xx)+FIN(x2).(xx-x1)
c  where x2-x1=1 on the FIN input grid.
c
c  JINF is the index of the first input point, or equivalently,
c  the number of points missing from the beginning of FIN,
c  similar to IFIRST in runlogs.
c
c  Typically we will have an infinite-resolution input spectrum FIN
c  calculated on a known grid DNU.  This will have been calculated
c  over a wider interval than the prescribed output window to allow
c  for the width of the ILS and for stretching.
c
c  Calculated Input spectrum is sampled at wavenumbers:
c     (JINF+J-1).DNU  J=1,NIN
c
c  And so covers points with wavenumbers
c     JINF.DNU  to  (JINF+NIN-1).DNU   (NIN points)
c
c  The output vector is obtained by re-sampling FIN
c  at the wavenumbers corresponding to the
c  already-read measured spectrum.
c     (KOUTF+K-1).RDEC.DNU   K=1,NOUT
c
c  Typically, we don't actually know the point spacing of the
c  measured spectrum -- it is determined from the fitting.
c  So RDEC might change slightly during iteration due to
c  frequency stretching after the measured spectrum has
c  been read. Since we don't want to re-read the measured
c  spectrum in again, slightly shifted in frequency, the
c  fitted window will deviate slightly from that specified
c  in cases with large frequency shifts/stretches.
c
c  The first output point must satisfy:
c      rdec*koutf >= jinf+(nhw-1+1/odec)
c  The last output point must satisfy:
c      rdec*(koutf+nout-1) <= jinf+nin-1-(nhw-1+1/odec)
c
c Implementation Notes:
c
c 1) Input vector (FIN) is assumed to be band-width limited
c   such that it can be perfectly interpolated/resampled
c   using a wide-enough SINC function.
c
c 2) All vectors are assumed to have equally spaced abscissae
c   that are integral multiples of their spacing and therefore
c   have a point (real or virtual) at the origin). The output
c   grid is larger than the input grid by a factor RDEC
c        FIN   has spacing of  1
c        OPER  has spacing of  1/ODEC  (ODEC is integer)
c        FOUT  has spacing of  RDEC    (RDEC is real)
c   The output vector can be under- or over-sampled wrt FIN.
c   RDEC can be <1 (interpolation) or >1 (decimation). In the
c   latter case, the output may not be band-width limited.
c
c 3) Since output grid is not typically an integer multiple
c    of the input grid, OPER is oversampled by an integer
c    factor ODEC with respect to FIN to allow more accurate
c    interpolation onto the non-integral output grid.
c
c 4) The length of the convolution operator OPER, is assumed to
c    be  NOP = 1+2*NHW*ODEC, which is always odd. OPER(1) and
c    OPER(NOP) are both zero, which allows all the dot products
c    to be performed with just 2*NHW operations, and allows use
c    of  OPER(1) or OPER(NOP) without risk of array-bound
c    violations of FIN.
c
c 5) In previous versions of this subroutine, benign array-bound
c    violations were allowed: e.g. when fr=0.0 the computation
c    of tt1 occasionally addressed an out-of-bound element of OPER,
c    or when fr=1.0 the computation of tt2 used an out-of-bound
c    element. These were benign because TT1 was multiplied by fr
c    and TT2 by (1-fr).
c        fout(kout)=fr*TT1+(1.0d0-fr)*TT2
c    For this reason, OPER was previously dimensioned (*) rather
c    than (1+2*NHW*ODEC). The latter would crash the program by
c    detecting and reporting array-bound violations. 
c    In the latest version, these violations are no longer allowed,
c    allowing OPER to be declared as OPER(1+2*NHW*ODEC). This
c    required a more complicated code for indexing the OPER array,
c    and also results in a slight reduction in the width of the
c    interpolable output vector.
c
c 6) The index of the first point of OPER, to be used in a dot
c    product, JOP+1, can only take on values from 1 to ODEC, which
c    means that the very last element of OPER is never used in the
c    evaluation of TT2. Only in the evaluation of TT1 can the first
c    element of OPER attain the value ODEC+1, which means that the
c    very last element of OPER will be used.
c        tt1=tt1+dprod(fin(jxx),oper(jop+2))
c        tt2=tt2+dprod(fin(jxx),oper(jop+1))
c
c 7) This subroutine is a replacement for newdec.f. Whereas newdec
c    used an initial fractional shift parameter to determine where
c    to start sampling of the input vector, regrid uses the integer
c    origin shifts of the input and output vectors.
c
c 8) When calling REGRID2, the user must specify the spectral
c    intervals (starting index, number of points) of both the
c    input and output vectors, together with the length of the
c    convolution operator. If these are inconsistent, in the
c    sense that the requested output vector cannot be fully
c    determined from the provided input vector and operator
c    without incurring an array bound violation, then REGRID
c    generates as much of the output vector as it can by
c    convolving with OPER, and for the remaining points it
c    performs a linear extrapolation from the left/right-most
c    two points of the input vector. When this happens, REGRID
c    returns non-zero values for the arguments NEXPL and NEXPR,
c    which signify the number of extrapolated output points to
c    the left and right of the input vector, respectively.
c
c    So output points NEXPL+1 to NOUT-NEXPR are computed
c    properly (by convolution with OPER) whereas points
c    1 to NEXPL and NOUT-NEXPR+1 to NOUT are computed
c    crudely (linear extrapolation).
c
c 9) Indices starting with J represent the input vector/spectrum
c    Indices starting with K represent the output vector/spectrum
c
c   Suppose we want
c   Fout (k+koutf) =  FIN[rdec*(j+jinf)]
c   regrid uses rdec & jinf to determine the sampling comb
c
c   Let sh=rdec*jinf-koutf
c     Fout(k) =  FIN[rdec*j+sh]
c   newdec used rdec & sh to determine the sampling comb.
c
c   Call requests output points from KOUTF to KOUTF+NOUT-1
c   Let KOUTL = KOUTF+NOUT-1
c   Points KOUTF to KOUTF-1+NEXPL are set to first good value (KOUTF+NEXPL).
c   Points KOUTL to KOUTL-1+NEXPL are set to last good value (KOUTL-NEXPR).

c--------------------------------------------------------------------------
      implicit none
      integer*4 nin,kout,nout,odec,nhw,jl,
     & jop,jinf,jinl,koutf,koutl,nexpl,nexpr,jxx

      real*4 fin(nin),
     & oper(1+2*nhw*odec),  ! causes benign array bound violations
     & fout(nout)

      real*8 
     & tt1,tt2,dvdot,dwid,rdec,xx,oxx,fr

c  Compute number of points that can't be computed by convolution
c  with OPER on the left and the right of the input vector.
      jinl=jinf+nin-1
      koutl=koutf+nout-1
      dwid=nhw-1+1/odec-1.0D-12
      if(rdec*koutf .lt. jinf+dwid) then
         nexpl=1+int((jinf+dwid)/rdec)-koutf
      else
         nexpl=0
      endif
c      write(*,*) 'jinf,dwid,koutf = ', jinf,dwid,koutf
c      write(*,*)rdec*koutf,' <',jinf+dwid,nexpl

      if(rdec*koutl .gt. jinl-dwid) then
         nexpr=koutl-int((jinl-dwid)/rdec)
      else
         nexpr=0
      endif

c      write(*,*)'Output units: ',koutf,(jinf+dwid)/rdec,nexpl
c      write(*,*)'Output units: ',koutl,(jinl-dwid)/rdec,nexpr

c  Convolution with OPER cannot be done for points KOUTF to
c  KOUTF-1+NEXPL.  So instead linearly extrapolate from two
c  nearest good values. This code is skipped if NEXPL=0
      do kout=1,min(nout,nexpl)
         xx=rdec*(koutf+kout-1)-jinf+1
         jl=max(int(xx),1)
         fr=xx-jl
         fout(kout)=sngl(dble(fin(jl))*(1-fr)+dble(fin(jl+1))*fr)
c         write(*,*)'l: kout,rdec*kout,xx,fout,fr = ',kout+koutf-1,
c     &   sngl(rdec*(kout+koutf-1)),sngl(xx-nhw+1),fout(kout),sngl(fr)
      end do

c  Main loop: over output points.
c  Compute function values (tt1, tt2) at sub-grid points that
c  bracket desired abscissa (xx), then linearly interpolate.
      do kout=1+nexpl,nout-nexpr
         xx=rdec*(koutf+kout-1)-jinf+1-nhw
         jxx=int(xx)+1
         oxx=odec*(jxx-xx)
         jop=int(oxx)
c         if(jop.ge.odec) then
c            jop=jop-odec
c            oxx=oxx-odec
c            jxx=jxx-1
c         endif
         fr=oxx-jop
         tt1=0.0d0
         if (fr.gt.0.0d0) tt1=dvdot(fin(jxx),1,oper(jop+2),odec,2*nhw)
         tt2=0.0d0
         if (fr.lt.1.0d0) tt2=dvdot(fin(jxx),1,oper(jop+1),odec,2*nhw)
c         do ii=1,2*nhw
c            tt1=tt1+dprod(fin(jxx),oper(jop+2))
c            tt2=tt2+dprod(fin(jxx),oper(jop+1))
c            jxx=jxx+1
c            jop=jop+odec
c         end do
         fout(kout)=sngl(fr*tt1+(1.0d0-fr)*tt2)  ! linear interpolation
c         write(*,*)'m: kout,rdec*kout,xx,fout,fr = ',kout+koutf-1,
c     &   sngl(rdec*(kout+koutf-1)),sngl(xx),fout(kout),sngl(fr)
      end do    !  do kout=1,nout

c  Convolution with OPER cannot be done for points KOUTL-1+NEXPL
c  to KOUTL.  So instead linearly extrapolate from nearest two good
c  values. This code is skipped if NEXPR=0
      do kout=max(1,nout-nexpr+1),nout
         xx=rdec*(koutf+kout-1)-jinf+1
         jl=min(int(xx),nin-1)
         fr=xx-jl
         fout(kout)=sngl(dble(fin(jl))*(1-fr)+dble(fin(jl+1))*fr)
c         write(*,*)'r: kout,rdec*kout,xx,fout,fr = ',kout+koutf-1,
c     &   sngl(rdec*(kout+koutf-1)),sngl(xx-nhw+1),fout(kout),sngl(fr)
      end do

      return
      end

      function dvdot(v1,inc1,v2,inc2,n)
      integer*4 i1,i2,j
      real*4 v1(*),v2(*)
      real*8 dvdot

      dvdot=0.0d0
      i1=1
      i2=1
      do j=1,n
         dvdot=dvdot+dprod(v1(i1),v2(i2))
         i1=i1+inc1
         i2=i2+inc2
      end do
      return
      end
