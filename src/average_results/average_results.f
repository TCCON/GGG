c  average_results.f
c  Averages results from different windows of the same gas 
c  after removing any window-related bias. Writes the
c  average values (column abundances) to the file runlog.xav
c  where "x" represents the geophysical quantity.
c
c  Arrays YOBS and YERR are physically 1-D for maximum flexibility,
c  but conceptually they are effectively 2-D arrays YOBS(nrow,nwin)
c  and are treated as such in the averaging_with XXX_bias subroutines
      implicit none
      include "../ggg_int_params.f"

      integer irow,jj,k,lnbc,mlabel,mval,navg,iav,lr,
     & lunr_xsw,lunw_xav,lunw_outlier,lunw_cew,lunw_rpt,
     & mwin,nwin,iwin,ngas,kgas,
     & mrow,
     & nauxcol,nlhead,ncol,jcol,icol,cwas,jav,
     & nrow,nss,loc_,loc,lwas
      parameter (lunr_xsw=51)      ! input file (.xsw)
      parameter (lunw_xav=54)      ! output file (.xav)
      parameter (lunw_cew=56)      ! wincomp file (.cew)
      parameter (lunw_rpt=57)      ! averages file (.rpt)
      parameter (lunw_outlier=58)  ! outlier file (.outlier)
      parameter (mwin=600)         ! Total number of columns/windows
      parameter (mrow=360000)      ! Max number of output records/spectra
      parameter (mval=20000000)    ! Max number of values (NROW * NCOL)
      parameter (mlabel=18000)     ! Max Number of column lable characters

      integer avindx(mgas+1), spectrum_flag,ntot(mgas)
      character
     & version*64,
     & gfit_version*64,
     & gsetup_version*64,
     & collate_version*64,
     & swfile*80,
     & avfile*80,
     & collabel*(mlabel),
     & sign(mrow)*1,
     & ftype*1,
     & clab(2*mwin+mauxcol)*24,
     & spectrum(mrow)*57, 
     & avlabel(mgas+1)*8,
     & data_fmt*40,
     & input_fmt*40,
     & output_fmt*40

      real*8 year(mrow),ty(mgas),t2(mgas),te(mgas),wt

      real*4
     & ymiss,rew(mrow),cew(mwin),tew,error_sigma,
     & yaux(mauxcol,mrow),
     & yobs(mval),yerr(mval),
     & ybar(mrow),eybar(mrow),
     & bias(mwin),ebias(mwin)

      version=
     &' average_results              Version 1.21    2012-10-05   GCT'
      write(*,*) version
      spectrum_flag=0
      loc=0

      if (iargc() == 0) then
         write(*,'(a)')
     &    'Enter name of .?sw file whose contents are to be averaged'
          read(*,'(a)') swfile
      elseif (iargc() == 1) then
          call getarg(1, swfile)
      else
          stop 'Usage: $gggpath/bin/average_results xswfile'
      endif
      lr=lnbc(swfile)
      avfile=swfile(:lr-2)//'av'
      ftype=swfile(lr-2:lr-2)

c  Read the entire contents of the .xsw disk file

      open(lunr_xsw,file=swfile,status='old')
c      open(lunw_xav,file=avfile,status='unknown')
      read(lunr_xsw,'(i2,i4,i7,i4)') nlhead,ncol,nrow,nauxcol
      if(nrow.gt.mrow) stop 'increase parameter mrow'
      read(lunr_xsw,'(a)') collate_version
      read(lunr_xsw,'(a)') gfit_version
      read(lunr_xsw,'(a)') gsetup_version
      read(lunr_xsw,'(8x,e12.5)') ymiss
      read(lunr_xsw,'(7x,a)') data_fmt
      read(lunr_xsw,'(a)') collabel
      if (index(collabel, 'spectrum') .gt. 0) spectrum_flag=1
      nwin=(ncol-nauxcol)/2

c      if (spectrum_flag .eq. 1) then
c         input_fmt='(a1,a38,1x,f13.8,NNf13.5,MMM(e12.4))'
c         write(input_fmt(18:19),'(i2.2)') nauxcol-1-spectrum_flag
c         write(input_fmt(26:28),'(i3.3)') ncol-nauxcol
c      else
c         input_fmt='(a1,f13.8,NNf13.5,MMM(e12.4))'
c         write(input_fmt(11:12),'(i2.2)') nauxcol-1-spectrum_flag
c         write(input_fmt(19:21),'(i3.3)') ncol-nauxcol
c      endif
      write(*,*) 'data fmt =  ',data_fmt
c      write(*,*) 'input format =  ',input_fmt
      input_fmt=data_fmt

      do irow=1,mrow
        if (spectrum_flag .eq. 1) then
           read(lunr_xsw,input_fmt,end=99)
     $     spectrum(irow),sign(irow),year(irow),
     $     (yaux(k,irow),k=3,nauxcol),
     $     (yobs(irow+nrow*(k-1)),yerr(irow+nrow*(k-1)),k=1,nwin)
        else
           read(lunr_xsw,input_fmt,end=99)
     $     sign(irow),year(irow),
     $     (yaux(k,irow),k=2,nauxcol),
     $     (yobs(irow+nrow*(k-1)),yerr(irow+nrow*(k-1)),k=1,nwin)
        endif
      end do  !  irow=1,mrow
99    close(lunr_xsw)
      if(nrow.ne.irow-1) stop 'NROW mismatch'
c
      open(lunw_cew,file=avfile(:lr)//'.cew',status='unknown')
      write(lunw_cew,'(2i3)') 2,4
      write(lunw_cew,'(a)')'  Window      Mean_Col   Std_Dev   Chi2/N'

      open(lunw_outlier,file=avfile(:lr)//'.outliers',status='unknown')
      call substr(collabel,clab,2*mwin+mauxcol,nss)
      if(nss.ne.ncol) stop 'NSS .NE. NCOL'
      write(*,*)'nrow=',nrow
      write(*,*)'ncol=',ncol
      write(*,*)'naux=',nauxcol
      write(*,*)'nwin=',nwin
      write(*,*)
c
c  Find location in collable of first _ (usually 'air_4444')
c  beyond the auxiliary variables.
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
      end do   !  do iwin=2,nwin
      ngas=kgas-1
      avindx(kgas)=iwin
      avlabel(kgas)=clab(icol)(:loc-1)

      write(*,*) '        Col1        Col2        Nwin   Gas'
      do kgas=1,ngas
         navg=avindx(kgas+1)-avindx(kgas)
         write(*,*)avindx(kgas),avindx(kgas+1)-1,navg,
     &   '   '//avlabel(kgas)
         jj=1+nrow*(avindx(kgas)-1)
         if(ftype.eq.'l' .or. ftype.eq.'v' .or. ftype.eq.'t') then
           call average_with_mul_bias(ymiss,nrow,navg,yobs(jj),yerr(jj),
     &     ybar,eybar,bias,ebias,rew,cew,tew)
           do irow=1,nrow     !   Report outliers
c              write(*,*)irow,ybar(irow),eybar(irow),rew(irow)
              jcol=avindx(kgas) 
              do jav=1,navg
                 jj=irow+nrow*(jcol-1)
                 if(yerr(jj).ne.ymiss) then
                 error_sigma=(yobs(jj)-ybar(irow)*bias(jav))/yerr(jj)
                 if(abs(error_sigma).gt.5.) 
     &           write(lunw_outlier,'(a12,f9.5,a22,i6,1x,a20,a11,a10)')
     &           ' Deviation =', error_sigma,' sigma for spectrum # ',
     &           nint(yaux(4+spectrum_flag,irow)),spectrum(irow)(:20),
     &           ' in window ',clab(nauxcol+2*jcol-1)
                 endif
                 jcol=jcol+1
              end do
           end do
         else
           call average_with_add_bias(ymiss,nrow,navg,yobs(jj),yerr(jj),
     &     ybar,eybar,bias,ebias,rew,cew,tew)
           do irow=1,nrow     !   Report outliers
              jcol=avindx(kgas) 
              do jav=1,navg
                jj=irow+nrow*(jcol-1)
                if(yerr(jj).ne.ymiss) then
                error_sigma=(yobs(jj)-ybar(irow)-bias(jav))/yerr(jj)
                if(abs(error_sigma).gt.5.) 
     &          write(lunw_outlier,'(a12,f9.5,a22,i6,2x,a11,a10)')
     &          ' Deviation =', error_sigma,' sigma for spectrum # ',
     &          nint(yaux(4+spectrum_flag,irow)),
     &          ' in window ',clab(nauxcol+2*jcol-1)
                endif
                jcol=jcol+1
             end do
         end do
         endif
c
c
c  Move the averages values back into the input array, overwriting
c  the earlier stuff that has already been averaged.
         call vmov( ybar,1,yobs(1+nrow*(kgas-1)),1,nrow)
         call vmov(eybar,1,yerr(1+nrow*(kgas-1)),1,nrow)
c
         if(navg.gt.1) then
         do iav=1,navg
            write(lunw_cew,'(a12,3f10.5)')
     &      clab(nauxcol+2*(avindx(kgas)+iav-1)-1),
     &      bias(iav),ebias(iav)*sqrt(float(nrow)),cew(iav)
         end do
         write(lunw_cew,*)
         endif

      end do   ! do kgas=1,ngas
      close(lunw_outlier)
      close(lunw_cew)

c
      if (spectrum_flag .eq. 1) then
         output_fmt='(a57,a1,f13.8,NNf13.5,MMM(1pe12.4))'
         write(output_fmt(15:16),'(i2.2)') nauxcol-1-spectrum_flag
         write(output_fmt(23:25),'(i3.3)') 2*ngas
      else
         output_fmt='(a1,f13.8,NNf13.5,MMM(1pe12.4))'
         write(output_fmt(11:12),'(i2.2)') nauxcol-1-spectrum_flag
         write(output_fmt(19:21),'(i3.3)') 2*ngas
      endif
c
c  Write averaged values to file
      open (lunw_xav,file=avfile,status='unknown')
      write(lunw_xav,'(i2,i4,i7,i4)')
     &      nlhead+1,nauxcol+2*ngas,nrow,nauxcol
      write(lunw_xav,'(a)') version(:lnbc(version))
      write(lunw_xav,'(a)') collate_version(:lnbc(collate_version))
      write(lunw_xav,'(a)') gfit_version(:lnbc(gfit_version))
      write(lunw_xav,'(a)') gsetup_version(:lnbc(gsetup_version))
      write(lunw_xav,'(a8,1pe12.4,a)') 'missing:',ymiss
      write(lunw_xav,'(a)') 'format='//output_fmt
      write(lunw_xav,'(a,4x,60(2x,a8,a14))') collabel(:loc_-4),
     & (avlabel(kgas),
     &  avlabel(kgas)(:lnbc(avlabel(kgas)))//'_error',kgas=1,ngas)
      do k=1,ngas
         ty(k)=0.0d0
         t2(k)=0.0d0
         te(k)=0.0d0
         ntot(k)=0
      end do
      do irow=1,nrow
         if (spectrum_flag .eq. 1) then   ! Spectrum names are included
            write(lunw_xav,output_fmt)
     &      spectrum(irow),sign(irow),year(irow),
     &      (yaux(k,irow),k=3,nauxcol),
     &      (yobs(irow+nrow*(k-1)),yerr(irow+nrow*(k-1)),k=1,ngas)
         else
c            write(lunw_xav,'(a1,f13.8,22f13.5,200(1pe12.4))')
            write(lunw_xav,output_fmt)
     &      sign(irow),year(irow),
     &      (yaux(k,irow),k=2,nauxcol),
     &      (yobs(irow+nrow*(k-1)),yerr(irow+nrow*(k-1)),k=1,ngas)
         endif
         do k=1,ngas
           if( yerr(irow+nrow*(k-1)).gt.1.0E-36) then
           wt=(1.0d0/yerr(irow+nrow*(k-1)))**2
           ty(k)=ty(k)+wt*yobs(irow+nrow*(k-1))
           t2(k)=t2(k)+wt*(dble(yobs(irow+nrow*(k-1))))**2
           te(k)=te(k)+wt
           ntot(k)=ntot(k)+1
           endif
c        write(*,*)irow,k,yobs(irow+nrow*(k-1)),yerr(irow+nrow*(k-1)),wt
         end do
      end do  ! do irow=1,nrow
      close(lunw_xav)
c
c  Write means for each gas and their standard errors and std deviations.
      open(lunw_rpt,file='average_results.rpt',status='unknown')
      write(lunw_rpt,*) '  k   ntot   gas         ybar       sterr
     &  stdev       Chi-2/N'
      do kgas=1,ngas
       write(lunw_rpt,'(i4,i8,2x,a,3(1pe12.4),0pf10.2)')kgas,ntot(kgas),
     &  avlabel(kgas), ty(kgas)/te(kgas),sqrt(ntot(kgas)/te(kgas)),
     & sqrt(t2(kgas)/te(kgas)-(ty(kgas)/te(kgas))**2),
     & sqrt(t2(kgas)/ntot(kgas)-ty(kgas)**2/te(kgas)/ntot(kgas))
      end do
      close(lunw_rpt)
      stop
      end
