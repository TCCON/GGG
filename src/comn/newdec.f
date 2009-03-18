      subroutine newdec(fin,nin,oper,nop,odec,rdec,sh,fout,nout)
c  Convolves input vector FIN(NIN) with the operator OPER(NOP). 
c  The result is placed in FOUT(NOUT). OPER is oversampled by an
c  integer factor ODEC with respect to FIN to allow interpolation.
c  The output vector can be under- or over-sampled wrt FIN by the
c  factor RDEC, which can be <1 (interpolation) or >1 (decimation).
c
c  INPUTS:
c      FIN(NIN)   R*4   discrete input function
c      NIN        I*4   number of elements of FIN that we are allowed to address
c      OPER(NOP)  R*4   operator to be convolved with FIN
c      NOP        I*4   number of elements in the operator OPER
c      ODEC       I*4   integer factor by which OPER is oversampled wrt FIN
c      RDEC       R*8   spacing of points in FOUT wrt FIN.
c      SH         R*8   initial fractional shift (in FIN grid points).
c      NOUT       I*4   Number of output points to be calculated
c
c  OUTPUTS:
c      FOUT(NOUT) R*4   Result of the convolution
c------------------------------------------------------------------------
c
c Implementation Notes:
c
c 1) All vectors are assumed to have equally spaced abscissae
c        FIN   has a spacing of  1
c        OPER  has a spacing of  1/ODEC
c        FOUT  has a spacing of  RDEC
c
c 2) The length of the convolution operator OPER, is assumed to equal
c    NOP = 1+2*NHW*ODEC, with OPER(1) and OPER(NOP) both zero.  This allows
c    all the dot products to be performed with just 2*NHW operations,
c    and allows us to start VDOT at OPER(KOP+0) or OPER(KOP+1)  without
c    risk of array-bound violations.
c
c 3) The index of the first point of OPER to be used in a dot product, KOP,
c    can only take on values from 1 to ODEC, which means that the very last
c    element of OPER is never used in the evaluation of TOT0.  Only in the
c    evaluation of TOT1 can the first element of OPER attain the value ODEC+1,
c    which means that the very last element of OPER will be used.
c
c 4) If OPER is a weakly apodized SINC function of period =< 1 and symmetrical
c    about its mid-point, OPER(1+ODEC*NHW), then NEWDEC will interpolate the
c    input vector, FIN, onto a new (and possibly asynchronous) grid such that 
c               FOUT(K) = FIN( (K-1)*RDEC+SH+1+NHW )
c    Note that FOUT is left-shifted with respect to FIN by NHW input points.
c    Also note that if OPER has a period >1, then NEWDEC will resample and
c    degrade/apodize FIN at the same time.
c
c 5) For a particular K, FIN will be addressed from INT( (K-1)*RDEC+SH+1+1 )
c                                               to  INT( (K-1)*RDEC+SH+1+2*NHW )
c    With K=1,NOUT, to avoid array-bound violations we must ensure that
c               INT( (K-1)*RDEC+SH+1+1 )             >= 1      for K=1
c          and  INT( (K-1)*RDEC+SH+1+2*NHW )         <= NIN    for K=NOUT
c    which reduce to
c               KLO = INT( SH+2 )                    >= 1      
c          and  KHI = INT( (NOUT-1)*RDEC+SH+1+2*NHW ) <= NIN 
c    Tests, performed at the beginning of NEWDEC, check these two conditions
c    and give warnings if they are not true. This is necessary because
c    although array-bound violations occurring within NEWDEC (including
c    the calls to VDOT) will result in run-time errors, array-bound violations
c    occurring within VDOT will not cause run-time errors because VDOT has
c    no information about the actual sizes of the arrays.
c--------------------------------------------------------------------------
      implicit none
      integer*4 kin,nin,kout,nout,nop,odec,kop,nhw,nn,klo,khi
      real*4 fin(nin),oper(nop),fout(nout),tot0,tot1,fr
      real*8 rdec,sh,xx

c      write(*,*)'NEWDEC:',nin,nop,odec,rdec,sh,nout
c
c Test for cases which will cause array-bound violations
      nhw=(nop-1)/odec/2
      if(nhw.lt.1) then
         write(*,*)' NEWDEC called with NOP < 1+2*ODEC',nop,odec,nhw
c         stop ' NEWDEC called with NOP < 1+2*ODEC'  ! Commented GCT 2009-03-18
         nhw=1
      endif
      klo=int( sh+2 )
      if(klo.lt.1) write(*,*)' NEWDEC warning: KLO < 1:',klo,1,sh
c      write(*,*)'nout,rdec,sh,nhw=', nout,rdec,sh,nhw
      khi=int( (nout-1)*rdec+sh+1+2*nhw )
      if(khi.gt.nin) write(*,*)' NEWDEC warning: KHI > NIN:',khi,nin
c      write(*,*) nout, rdec, sh, nhw,khi,nin
c
c  Main loop
      xx=odec*(sh+1)
      do kout=1,nout
c
c  Interpolate linearly between two nearest operator points
         nn=int(xx)
         fr=xx-nn
         kin=1+nn/odec
         kop=odec-mod(nn,odec)
         if(fr.gt.0.0) call vdot(fin(kin),1,oper(kop+0),odec,tot0,2*nhw)
         if(fr.lt.1.0) call vdot(fin(kin),1,oper(kop+1),odec,tot1,2*nhw)
         fout(kout)=(1.0-fr)*tot1+fr*tot0
c
c  Use nearest point (faster, but prone to cause oscillation)
c         nn=nint(xx)-1
c         kin=1+nn/odec
c         kop=odec-mod(nn,odec)
c         call vdot(fin(kin),1,oper(kop+0),odec,fout(kout),2*nhw)
c
         xx=xx+rdec*odec
      end do    !  do kout=1,nout
      return
      end
