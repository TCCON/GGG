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
      integer irow,jj,k,lnbc,mlabel,mval,navg,iav,
     & lr,lunr,lunw,luno,lunt,mwin,nwin,iwin,ngas,kgas,
     & mrow,mauxcol,nauxcol,nlhead,mgas,ncol,jcol,icol,cwas,jav,
     & nrow,nss,locnaux,loc,lwas
      parameter (lunr=12)      ! input file (.xsw)
      parameter (lunw=14)      ! output file (.xav)
      parameter (lunt=15)      ! wincomp file (.xav)
      parameter (luno=16)      ! outlier file (.xsw)
      parameter (mgas=80)      ! Max number of gases
      parameter (mwin=600)     ! Total number of columns/windows
      parameter (mrow=240000)  ! Max number of output records/spectra
      parameter (mval=12000000) ! Max number of values (NROW * NCOL)
      parameter (mauxcol=25)   ! Number of auxiliary parameters/columns
      parameter (mlabel=18000) ! Max Number of column lable characters

      integer avindx(mgas+1), naux, nchar
      character
     & gfit_version*64,gsetup_version*64,
     & collabel*(mlabel),swfile*80,avfile*80,
     & sign(mrow)*1,ftype*1,
     & clab(2*mwin+mauxcol)*17,
     & collate_version*64,ar_version*64,
     & spectrum(mrow)*38, 
     & avlabel(mgas+1)*8,
     & input_fmt*40, output_fmt*40

      real*8 year(mrow)

      real*4
     & ymiss,rew(mrow),cew(mwin),tew,error_sigma,
     & yaux(mauxcol,mrow),
     & yobs(mval),yerr(mval),
     & ybar(mrow),eybar(mrow),
     & bias(mwin),ebias(mwin)

      ar_version=
     &' average_results              Version 1.1.2   2009-11-07   GCT'
      write(*,*) ar_version
      nchar=0
      loc=0

      write(*,'(a)')
     & 'Enter name of .?sw file whose contents are to be averaged'
      read(*,'(a)')swfile
      lr=lnbc(swfile)
      avfile=swfile(:lr-2)//'av'
      ftype=swfile(lr-2:lr-2)

c  Read the entire contents of the .xsw disk file

      open(lunr,file=swfile,status='old')
      open(lunw,file=avfile,status='unknown')
      read(lunr,'(i2,i4,i7,i4)') nlhead,ncol,nrow,nauxcol
      if(nrow.gt.mrow) stop 'increase parameter mrow'
      read(lunr,'(a)') collate_version
      read(lunr,'(a)') gfit_version
      read(lunr,'(a)') gsetup_version
      read(lunr,*) ymiss
      read(lunr,'(a)') collabel
      if (index(collabel, 'Spectrum') .gt. 0) nchar=1
      naux=nauxcol+nchar         ! ncol includes spectrum name
      nwin=(ncol-naux)/2

      if (nchar .eq. 1) then
         input_fmt='(a1,(a38,1x),f13.8,NNf13.5,800(e12.4))'
         write(input_fmt(20:21),'(i2.2)') nauxcol-1
      else
         input_fmt='(a1,f13.8,NNf13.5,800(e12.4))'
         write(input_fmt(11:12),'(i2.2)') nauxcol-1
      endif

      do irow=1,mrow
        if (nchar .eq. 1) then
           read(lunr,input_fmt,end=99)
     $     sign(irow),spectrum(irow),year(irow),
     $     (yaux(k,irow),k=2,nauxcol),
     $     (yobs(irow+nrow*(k-1)),yerr(irow+nrow*(k-1)),k=1,nwin)
        else
           read(lunr,input_fmt,end=99)
     $     sign(irow),year(irow),
     $     (yaux(k,irow),k=2,nauxcol),
     $     (yobs(irow+nrow*(k-1)),yerr(irow+nrow*(k-1)),k=1,nwin)
        endif
      end do  !  irow=1,mrow
99    close(lunr)
      if(nrow.ne.irow-1) stop 'NROW mismatch'
c
      open(lunt,file=avfile(:lr)//'.cew',status='unknown')
      write(lunt,'(2i3)') 2,4
      write(lunt,'(a)')' Window     Mean_Col   Std_Err   Chi2/N'

      open(luno,file=avfile(:lr)//'.outliers',status='unknown')
      call substr(collabel,clab,2*mwin+mauxcol,nss)
      if(nss.ne.ncol) stop 'NSS .NE. NCOL'
      write(*,*)nrow,nwin,nss
      locnaux=index(clab(naux+1),'_')
      locnaux=index(collabel,clab(naux+1)(:locnaux-1))-1
      cwas=naux-1
      icol=naux+1
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

      do kgas=1,ngas
         navg=avindx(kgas+1)-avindx(kgas)
         write(*,*)avindx(kgas),avindx(kgas+1)-1,navg,' '//avlabel(kgas)
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
     &          write(luno,'(a12,f9.5,a22,i6,a11,a10)')
     &          ' Deviation =', error_sigma,' sigma for spectrum # ',
     &          nint(yaux(4,irow)),' in window ',clab(naux+2*jcol-1)
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
     &          write(luno,'(a12,f9.5,a22,i6,a11,a10)')
     &          ' Deviation =', error_sigma,' sigma for spectrum # ',
     &          nint(yaux(4,irow)),' in window ',clab(naux+2*jcol-1)
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
            write(lunt,'(a10,3f10.5)')
     &      clab(naux+2*(avindx(kgas)+iav-1)-1),
     &      bias(iav),ebias(iav),cew(iav)
         end do
         write(lunt,*)
         endif

      end do   ! do kgas=1,ngas

c

      if (nchar .eq. 1) then
         output_fmt='(a1,(a38,1x),f13.8,NNf13.5,200(1pe12.4))'
         write(output_fmt(20:21),'(i2.2)') nauxcol-1
      else
         output_fmt='(a1,f13.8,NNf13.5,200(1pe12.4))'
         write(output_fmt(11:12),'(i2.2)') nauxcol-1
      endif
c
c  Write averaged values to file
      open (lunw,file=avfile,status='unknown')
      write(lunw,'(i2,i4,i7,i4)') nlhead+1,naux+2*ngas,nrow,nauxcol
      write(lunw,'(a)') ar_version(:lnbc(ar_version))
      write(lunw,'(a)') collate_version(:lnbc(collate_version))
      write(lunw,'(a)') gfit_version(:lnbc(gfit_version))
      write(lunw,'(a)') gsetup_version(:lnbc(gsetup_version))
      write(lunw,'(1pe12.4,a)') ymiss,'   ! missing value'
      write(lunw,'(a,60(a8,2x,a14))') collabel(:locnaux),
     & (avlabel(kgas)(:lnbc(avlabel(kgas))),
     &  avlabel(kgas)(:lnbc(avlabel(kgas)))//'_error',kgas=1,ngas)
      do irow=1,nrow
        if (nchar .eq. 1) then
           write(lunw,output_fmt)
     &     sign(irow),spectrum(irow),year(irow),
     &     (yaux(k,irow),k=2,nauxcol),
     &     (yobs(irow+nrow*(k-1)),yerr(irow+nrow*(k-1)),k=1,ngas)
        else
c           write(lunw,'(a1,f13.8,22f13.5,200(1pe12.4))')
           write(lunw,output_fmt)
     &     sign(irow),year(irow),
     &     (yaux(k,irow),k=2,nauxcol),
     &     (yobs(irow+nrow*(k-1)),yerr(irow+nrow*(k-1)),k=1,ngas)
        endif
      end do
      close(lunw)
      close(luno)
      stop
      end
