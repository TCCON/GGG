c  average_results.f
c  Reads the runlog.xsw input file and averages results
c  from consecutive windows of the same gas after removing
c  any window-related bias. Writes the average values
c  (e.g. column abundances) to the file runlog.xav where
c  "x" represents the geophysical quantity:
c     v = vertical column
c     l = line of sight column
c     t = VMR Scale Factor
c
c  Input file
c     runlog.xsw or runlog.xsw.ada
c
c  Output files:
c     runlog.xav
c     runlog.xav.cew
c     runlog.xsw.transpose
c     runlog.xav.outliers
c     runlog.xav.rew
c     runlog.xav.cc
c     average_results.rpt
c
c  Arrays YOBS and YERR are physically 1-D for flexibility,
c  but conceptually they are 2-D arrays YOBS(nrow,nwin) and
c  are treated 2-D in the averaging_with_xxx_bias subroutines.
c
c  The program has two modes of operation. If the 6'th row of
c  the .xsw input file header begins with "sf=", then the
c  subsequent biases values are assumed in the correction,
c  This is achieved by setting mit=1, suppressing iteration.
c  If no sf= line is found, then mit is set to 25 and the
c  program iterates to find its own bias values, as it did
c  previously. Thus the program can be made to function as
c  previously simply by removing the sf= line from the header
c  (and reducing nlhead from 8 to 7). In other words, the
c  program will process old-syle .xsw files in the old way
c  (iterating) and new-syle .xsw files by applying the
c  sf-values given in the header.

      implicit none
      include "../gfit/ggg_int_params.f"
      include "../comn/postproc_params.f"

      integer irow,j,jj,k,kk,kflag,lnbc,mlabel,mval,mavg,navg,kav,jav,
     & lr,lunr_xsw,lunw_xav,lunw_outlier,eflag,
     & lunw_cew,lunw_rew,lunw_cc,lunw_rpt,
     & mwin,nwin,iwin,ngas,kgas,mrow,idum,
     & nauxcol,nlhead,ncol,kcol,jcol,icol,cwas,
     & nrow,nss,loc_,loc,lwas,mit,naddn,ihead, kxo2
      logical isada, force_iter
      logical isclose_s            ! function that tests whether real*4
                                   ! values are within floating point error 
                                   ! of each other

      parameter (lunr_xsw=51)      ! input file (.xsw)
      parameter (lunw_xav=54)      ! output file (.xav)
      parameter (lunw_cew=55)      ! wincomp file (.cew)
      parameter (lunw_rew=56)      ! .rew file (.cew)
      parameter (lunw_cc=57)       ! .rew file (.cew)
      parameter (lunw_rpt=58)      ! averages file (.rpt)
      parameter (lunw_outlier=59)  ! outlier file (.outlier)
      parameter (mwin=621)         ! Total number of columns/windows
      parameter (mavg=102)         ! Max # of columns/windows to be averaged
      parameter (mrow=360000)      ! Max number of output records/spectra
      parameter (mval=30000000)    ! Max number of values (NROW * NCOL)
      parameter (mlabel=18000)     ! Max Number of column lable characters

      integer avindx(mgas+1), spectrum_flag,nused(mgas)
      character
     & version*64,
     & sfstring*2048,
     & swfile*80,
     & avfile*80,
     & collabel*(mlabel),
     & col1(mrow)*1,
     & ftype*1,
     & clab(2*mwin+mauxcol)*24,
     & spectrum(mrow)*57, 
     & avlabel(mgas+1)*8,
     & data_fmt*40,
     & input_fmt*40,
     & output_fmt*40,
     & addn_lines(maddln)*(mcharhead),
     & fiter_cl*6

      real*8 year(mrow),ty(mgas),t2(mgas),te(mgas),wt,
     & dtiny,covmat(mavg*mavg),ymiss_dbl,tt
      real*8 wlimit

      real*4
     & yobs(mval),yerr(mval), ! Data from inout file
     & ccbar(mavg),yy,small,
     & ymiss,rew(mrow,mwin),cew(mwin),tew,error_sigma,
     & yaux(mauxcol,mrow),
     & ybar(mrow),eybar(mrow),
     & bias(mwin),ebias(mwin)

      idum=mcolvsw    ! prevent compiler warning (unused parameter)
      idum=mcolvav    ! prevent compiler warning (unused parameter)
      idum=mfilepath  ! prevent compiler warning (unused parameter)
      idum=mlev       ! prevent compiler warning (unused parameter)
      idum=mrow_qc    ! prevent compiler warning (unused parameter)
      idum=mspeci     ! prevent compiler warning (unused parameter)
      idum=mvmode     ! prevent compiler warning (unused parameter)
      idum=ncell      ! prevent compiler warning (unused parameter)
      idum=nchar      ! prevent compiler warning (unused parameter)

      version=
     &' average_results          Version 1.37    2020-07-31   GCT,JLL'
      write(*,*) version
      spectrum_flag=0
      loc=0
      dtiny=1.D-88
      small=1.0E-18

      if (iargc() == 0) then
         write(*,'(a)')
     &   'Enter name of .?sw file whose contents are to be averaged'
         read(*,'(a)') swfile
      elseif (iargc() == 1) then
         call getarg(1, swfile)
      elseif (iargc() == 2) then
         call getarg(1, swfile)
         call getarg(2, fiter_cl)
      else
         stop 'Usage: $gggpath/bin/average_results xswfile'
      endif

      force_iter = fiter_cl .eq. '--iter'

      lr=lnbc(swfile)
c JLL: average_results may be getting an .xsw or .xsw.ada file
c it doesn't matter which except that we need to change the right
c part of the file name (sw -> av)
      isada = .false.
      if (swfile(lr-2:lr) .eq. 'ada') then
c       write(*,*) 'ada file'
        write(*,*) swfile(:lr)
        isada = .true.
        avfile=swfile
        avfile(lr-5:lr-4) = 'av'
        ftype=swfile(lr-6:lr-6)
      else
c       write(*,*) 'non-ada file'
        write(*,*) swfile(:lr)
        avfile=swfile(:lr-2)//'av'
        ftype=swfile(lr-2:lr-2)
      end if

c Read the xsw file, starting with the header
      open(lunr_xsw,file=swfile,status='old')
      call read_postproc_header(lunr_xsw, nlhead, ncol, nrow, nauxcol,
     & ymiss_dbl, data_fmt, addn_lines, naddn)

      ymiss = sngl(ymiss_dbl)      
      if(nrow.gt.mrow) stop 'increase parameter mrow'
      nwin = (ncol-nauxcol)/2

c go ahead and insert this program's version into the additional lines
      do ihead=naddn,1,-1
        addn_lines(ihead+1) = addn_lines(ihead)
      end do
      addn_lines(1) = version
      naddn = naddn + 1

c Initially, assume that we will need to iteratively find the window-to
c window biases. Set the number of iterations to 25 and initialize these
c biases to 1.
      mit = 25
      do jcol=1,nwin
        bias(jcol) = 1.0
      end do

c Search the header for a line containing "sf=", if that is
c present, then we use the scale factors defined in there.
c Once done, remove the scale factors from the header - they
c don't make sense post-averaging because they won't line up
c with the columns.
      do ihead=1,naddn
        if (index(addn_lines(ihead),'sf=') .gt. 0 
     &      .and. .not. force_iter) then
          mit = 1
          read(addn_lines(ihead), '(3x,400f8.3)') 
     & (bias(jcol),jcol=1,nwin)

          do j=ihead+1,naddn
            addn_lines(j-1) = addn_lines(j)
          end do
          naddn = naddn -1
          goto 30
        endif
      end do
 30   continue
c     write(*,*) 'mit = ', mit

      read(lunr_xsw,'(a)') collabel
c      write(*,*) collabel(:256)
      if (index(collabel, 'spectrum') .gt. 0) spectrum_flag=1
c      write(*,'(a,400f8.3)')'bias=',(bias(jcol),jcol=1,(ncol-nauxcol)/2)
      input_fmt=data_fmt

c  Read the entire contents of the .xsw disk file

      do irow=1,mrow
         if (spectrum_flag .eq. 1) then
            read(lunr_xsw,input_fmt,end=99) spectrum(irow),
     &      col1(irow),year(irow),(yaux(k,irow),k=3,nauxcol),
     $      (yobs(irow+nrow*(k-1)),yerr(irow+nrow*(k-1)),k=1,nwin)
         else
            read(lunr_xsw,input_fmt,end=99)
     $      col1(irow),year(irow),(yaux(k,irow),k=2,nauxcol),
     $      (yobs(irow+nrow*(k-1)),yerr(irow+nrow*(k-1)),k=1,nwin)
         endif

cc   Prevent the spectra unused by diurnret from influencing averaging
c         if(col1(irow).eq.'-') then
c            do k=1,nwin
c               yerr(irow+nrow*(k-1))=ymiss
c            end do
c         endif

      end do  !  irow=1,mrow
99    close(lunr_xsw)
      if(nrow.ne.irow-1) stop 'NROW mismatch'
c
      open(lunw_cew,file=avfile(:lr)//'.cew',status='unknown')
      write(lunw_cew,'(2i3)') 2,5
      write(lunw_cew,'(a)')
     &'  Window           VSF    VSF_error   Chi2/N    CC_bar'

c  Pre-delete existing outlier file in case no outliers are found.
c  (it won't be over-written). Then open new empty outlier file.
      open(lunw_outlier,file=avfile(:lr)//'.outliers',status='unknown')
      close(lunw_outlier,status='delete')
      open(lunw_outlier,file=avfile(:lr)//'.outliers',status='new')

c  Decompose collabel string into vector of substrings (clab)
      call substr(collabel,clab,2*mwin+mauxcol,nss)
      if(nss.ne.ncol) then
         write(*,*) 'NSS,NCOL=',nss,ncol
         stop 'NSS .NE. NCOL'
      endif
c      write(*,*)'nss =',nss
c      write(*,*)'nrow=',nrow
c      write(*,*)'ncol=',ncol
c      write(*,*)'naux=',nauxcol
c      write(*,*)'nwin=',nwin
c      write(*,*)

c JLL 2020-06-04: we'll need the xo2 error to add in if operating on 
c a .vsw.ada file. We'll just find it always b/c this should be fast
      kxo2 = 0
      do iwin=nauxcol+2,2*nwin+nauxcol,2
c        write(*,*) iwin, clab(iwin), index(clab(iwin),'xo2'),
c     & index(clab(iwin), 'error')
        if (index(clab(iwin), 'xo2').gt.0) kxo2 = (iwin-nauxcol)/2
      end do

      if (kxo2 .eq. 0 .and. isada) 
     & stop 'xo2_error column not found'
c
c  Find location in collabel of first underscore
c  (usually in 'luft_nnnn') beyond the auxiliary variables.
      loc_=index(clab(nauxcol+1),'_')
      loc_=index(collabel,clab(nauxcol+1)(:loc_-1))
      cwas=nauxcol-1
      icol=nauxcol+1
      lwas=1
      kgas=1
      do iwin=1,nwin
         loc=index(clab(icol),'_')
         if(loc.eq.0) loc=index(clab(icol),'-')
         if(clab(cwas)(:lwas-1).ne.clab(icol)(:loc-1)) then  ! new gas
            avindx(kgas)=iwin
            avlabel(kgas)=clab(icol)(:loc-1)
            kgas=kgas+1  ! index of window/column to write averages
            if(kgas.gt.mgas) stop 'kgas.gt.mgas'
         endif
         cwas=icol
         lwas=loc
         icol=icol+2
      end do   !  do iwin=1,nwin
      ngas=kgas-1
      avindx(kgas)=iwin
      avlabel(kgas)=clab(icol)(:loc-1)

      open(lunw_cc,file=avfile(:lr)//'.cc',status='unknown')
      open(lunw_rew,file=avfile(:lr)//'.rew',status='unknown')
      write(lunw_rew,'(2i3)') 4,ngas+1
      write(lunw_rew,'(a)') ' Row Error Weights (.rew values).'
      write(lunw_rew,'(a)') ' Gases with only one window have zeros.'
      write(lunw_rew,'(a9,48a13)')' ispec   ',
     & (avlabel(kgas),kgas=1,ngas)

c      write(*,*) '      kgas   indx(k)  indx(k+1)-1   navg   Gas'
      do kgas=1,ngas
         navg=avindx(kgas+1)-avindx(kgas)
c         write(*,'(4i10,a)')kgas,avindx(kgas),avindx(kgas+1)-1,navg,
c     &   '   '//avlabel(kgas)
         jj=1+nrow*(avindx(kgas)-1) ! pointer to start of kgas data

c  Write out lower triangle of correlation matrix
         if(navg.le.mavg) then
c           write(*,*) 'Calling covariance_matrix',navg,nrow,jj
            call covariance_matrix
     &       (ymiss,2,navg,nrow,yobs(jj),yerr(jj),covmat)
            write(lunw_cc,*)
            write(lunw_cc,*)'Correlation Coefficients:'
            write(lunw_cc,*)avindx(kgas),avindx(kgas+1)-1,navg,
     &      '   '//avlabel(kgas)
            write(lunw_cc,'(1x,20i7)') (jav,jav=1,navg)
            do kav=1,navg
               tt=0.0d0
               do jav=1,navg
                  tt=tt+covmat(kav+navg*(jav-1))/
     &            dsqrt(dtiny+covmat(kav+navg*(kav-1))*
     &            covmat(jav+navg*(jav-1)))
               end do
               ccbar(kav)=sngl(tt/navg)
               write(lunw_cc,'(i3,103f7.3)') kav,
     &         (covmat(kav+navg*(jav-1))/
     &         dsqrt(dtiny+covmat(kav+navg*(kav-1))*
     &         covmat(jav+navg*(jav-1))),jav=1,navg),ccbar(kav)
            end do
         endif   !  if(navg.le.mavg) then

c  Find ybar of the navg observations and the bias of each window
         if(ftype.eq.'l' .or. ftype.eq.'v' .or. ftype.eq.'t') then
c   error_sigma=(yobs(jj)-ybar(irow)*bias(jav))/yerr(jj)
            do jav=1,navg
c               bias(avindx(kgas)-1+jav)=1.0
               ebias(avindx(kgas)-1+jav)=5.
            end do
            call average_with_mul_bias(mit,ymiss,nrow,navg,yobs(jj),
     &      yerr(jj),ybar,eybar,bias(avindx(kgas)),ebias(avindx(kgas)),
     &      rew(1,kgas),cew,tew)
            kflag=1
         else
c   error_sigma=(yobs(jj)-ybar(irow)-bias(jav))/yerr(jj)
            do jav=1,navg
               bias(avindx(kgas)-1+jav)=0.0
               ebias(avindx(kgas)-1+jav)=5.
            end do
            call average_with_add_bias(mit,ymiss,nrow,navg,yobs(jj),
     &      yerr(jj),ybar,eybar,bias(avindx(kgas)),ebias(avindx(kgas)),
     &      rew(1,kgas),cew,tew)
            kflag=0
         endif  ! if(ftype.eq.'l' .or. ftype.eq.'v' .or. ftype.eq.'t')

         if(navg.gt.1) then
            do jav=1,navg
               write(lunw_cew,'(a14,4(1x,f9.4))') 
     &         clab(nauxcol+2*(avindx(kgas)+jav-1)-1),
     &         bias(avindx(kgas)-1+jav),ebias(avindx(kgas)-1+jav)*
     &         sqrt(float(nrow)),wlimit(dble(cew(jav)), 'f9.4'),
     &         ccbar(jav)
            end do
            write(lunw_cew,*)
         endif

         if(isada .and. avindx(kgas).ne.kxo2) then
c  JLL 2020-06-04: We've modified the error propagation now that
c  individual windows are airmass corrected. When averaging xgas 
c  values, we don't want the O2 error included in the xgas_window
c  error for two reasons:
c       1. If the O2 error is large relative to the gas error it will
c          dominate the window-to-window error weights, which is 
c          undesireable.
c       2. With >1 window, the O2 error will be correlated between
c          the windows (so it would be a systematic error) but 
c          when calculating the average's error, it is treated as
c          random, which is incorrect.
c  The solution is to add in the O2 error here, which avoids both
c  problems. Specifically we add X_i * xo2_error / 0.2095 in quadrature
c  to the current gas average error. The factor of 0.2095 converts the
c  xo2 error from a fraction of the air column to a fraction of the
c  O2 column.
            do irow=1,nrow
c               Pointer to the corresponding xo2 error value
                kk = nrow*(kxo2-1) + irow

c               Only if all quantities needed for this calculation 
c               are not missing values do we compute it. If ANY
c               are missing, then we set the average error to a
c               missing value. In most cases, it will already be,
c               but may as well be sure.
c                write(*,*) avlabel(kgas), 'irow=', irow, 'ymiss=', 
c     & ymiss, 'ybar=', ybar(irow), 'eybar=', eybar(irow), 'yerr=', 
c     & yerr(kk), 'tew=', tew, 'rew(1)=', rew(1,kgas)
                if( isclose_s(eybar(irow), ymiss)
     &              .or. isclose_s(ybar(irow), ymiss)
     &              .or. isclose_s(yerr(kk), ymiss) ) then

                    eybar(irow) = ymiss

                else
                    eybar(irow) = sqrt(eybar(irow)**2 +
     & (ybar(irow)*yerr(kk)/0.2095)**2)
                endif
c                write(*,*) 'eybar after=', eybar(irow)
            end do ! do irow=1,nrow
         else
c            write(*,*) 'Not adding O2 error for ', avlabel(kgas)
         endif  ! if(isada)

         eflag=0
         do irow=1,nrow     !   Report outliers
            if(abs(eybar(irow)).gt.ymiss) eflag=eflag+1
c            jcol=avindx(kgas) 
            if(navg.gt.1) then
               do jav=1,navg
                  jj=irow+nrow*((avindx(kgas)+jav-1)-1)
                  if(.not. isclose_s(yerr(jj), ymiss)) then
                     yy=kflag*ybar(irow)*bias(avindx(kgas)-1+jav)-
     &               (kflag-1)*(ybar(irow)+bias(avindx(kgas)-1+jav))
                     error_sigma=(yobs(jj)-yy)/yerr(jj)
                     if(abs(error_sigma).gt.6.) write(lunw_outlier,
     &               '(a12,f9.5,a22,i6,a11,a10,2f9.5)')' Deviation =',
     &               error_sigma,' sigma for spectrum # ',
     &               nint(yaux(4+spectrum_flag,irow)),' in window ',
     &               clab(nauxcol+2*(avindx(kgas)+jav-1)-1),
     &               yobs(jj),yerr(jj)
                  endif
c                  jcol=jcol+1
               end do
            endif
         end do     !  do irow=1,nrow     !   Report outliers
         if(eflag.gt.0) write(*,*) 'Warning: eflag = ',eflag
c
c  Move the average values values back into the input array,
c  overwriting the earlier stuff that has already been averaged.
c  This stores the average values prior to them all being written
c  to the .xav file at the end, but over-writes yobs and yerr.
         call vmov( ybar,1,yobs(1+nrow*(kgas-1)),1,nrow)
         call vmov(eybar,1,yerr(1+nrow*(kgas-1)),1,nrow)
      end do   ! do kgas=1,ngas
      close(lunw_outlier)
      close(lunw_cew)

c
      if (spectrum_flag .eq. 1) then
         output_fmt='(a57,a1,f13.8,NNf13.5,MMM(1pe13.5))'
         write(output_fmt(15:16),'(i2.2)') nauxcol-1-spectrum_flag
         write(output_fmt(23:25),'(i3.3)') 2*ngas
      else
         output_fmt='(a1,f13.8,NNf13.5,MMM(1pe13.5))'
         write(output_fmt(11:12),'(i2.2)') nauxcol-1-spectrum_flag
         write(output_fmt(19:21),'(i3.3)') 2*ngas
      endif
c
c  Write averaged values to file
      open (lunw_xav,file=avfile,status='unknown')
      call write_postproc_header(lunw_xav, nauxcol+2*ngas, nrow, 
     & nauxcol, dble(ymiss), output_fmt, addn_lines, naddn, 0)
      write(lunw_xav,'(a,4x,80(2x,a8,a14))') collabel(:loc_-4),
     & (avlabel(kgas),
     &  avlabel(kgas)(:lnbc(avlabel(kgas)))//'_error',kgas=1,ngas)
      do k=1,ngas
         ty(k)=0.0d0
         t2(k)=0.0d0
         te(k)=0.0d0
         nused(k)=0
      end do
      do irow=1,nrow
         write(lunw_rew,'(i9,444(1pe13.6))') irow,
     &   (rew(irow,jcol),jcol=1,ngas)
         if (spectrum_flag .eq. 1) then   ! Spectrum names are included
            write(lunw_xav,output_fmt)
     &      spectrum(irow),col1(irow),year(irow),
     &      (yaux(k,irow),k=3,nauxcol),
     &      (yobs(irow+nrow*(k-1)),yerr(irow+nrow*(k-1)),k=1,ngas)
         else
c            write(lunw_xav,'(a1,f13.8,22f13.5,200(1pe12.4))')
            write(lunw_xav,output_fmt)
     &      col1(irow),year(irow),
     &      (yaux(k,irow),k=2,nauxcol),
     &      (yobs(irow+nrow*(k-1)),yerr(irow+nrow*(k-1)),k=1,ngas)
         endif
         do k=1,ngas
            if( yerr(irow+nrow*(k-1)).gt.1.0E-36) then
               wt=(1.0d0/yerr(irow+nrow*(k-1)))**2
               ty(k)=ty(k)+wt*yobs(irow+nrow*(k-1))
               t2(k)=t2(k)+wt*(dble(yobs(irow+nrow*(k-1))))**2
               te(k)=te(k)+wt
               nused(k)=nused(k)+1
            endif
c            write(*,*)irow,k,yobs(irow+nrow*(k-1)),yerr(irow+nrow*(k-1)),wt
         end do  ! k=1,ngas
      end do  ! do irow=1,nrow
      close(lunw_xav)
      close(lunw_rew)
      close(lunw_cc)
c
c  Write means for each gas and their standard errors and std deviations.
      open(lunw_rpt,file='average_results.rpt',status='unknown')
      write(lunw_rpt,*)   '  k   nused   nrow    gas         ybar    
     & sterr      stdev      Chi-2/N'
      do kgas=1,ngas
         write(lunw_rpt,'(i4,2i8,2x,a,3(1pe12.4),0pf10.2)')kgas,
     &   nused(kgas),nrow,avlabel(kgas),ty(kgas)/te(kgas),
     &   sqrt(nused(kgas)/te(kgas)),
     &   sqrt(t2(kgas)/te(kgas)-(ty(kgas)/te(kgas))**2),
     &   dsqrt(t2(kgas)/nused(kgas)-ty(kgas)**2/te(kgas)/nused(kgas))
      end do
      close(lunw_rpt)
      stop
      end
