c  Program: apply_ghost_correction.f
c
c  Purpose: To apply the (airmass-independent) ghost correction to the data
c  in the selected runlog.vav.ada.aia file
c
c  yout(i,j,t) = GCF(j,t) + yin(i,j,t)
c
c  where
c   yin(i,j) is the DMF of the j'th gas from the i'th spectrum
c   yout(i,j) is the DMF of the j'th gas from the i'th spectrum
c   GCF(j) is the Ghost Correction Factor for the j'th gas at time t
c   which is read from the file ~/ggg/tccon/bi_ghost_corr.dat
c
c  Input Files:
c       runlog.vav.ada.aia 
c       sitename_ghost_corr.dat
c
c  Output Files:
c       runlog.gaa.vav.ada.aia
c       
      implicit none
      include "../gfit/ggg_int_params.f"
      include "../comn/postproc_params.f"

      integer*4 lunr,luns,lunw,ncoml,ncol,mcol,kcol,icol,j,
     & kgas,lnbc,irow,naux,ngas,nrow,li,k,ngrow,ngas1,kk,mrow
      parameter (lunr=14,luns=15,lunw=16,mcol=150,mrow=5)
      character header*800,headarr(mcol)*20,gasname(mgas)*20,
     & gggdir*(mpath),inputfile*40,outputfile*40, version*63,gaserr*32,
     & filename*(mpath+40)
      real*8 yrow(mcol),gcf(mrow,mgas),gcfe(mrow,mgas),scl(mrow,mgas),
     & cf(mcol,mgas)
      logical fileexists

      integer*4 endyear(mgas),endday(mgas)

      character outfmt*42,specname*(nchar),site*2,dl*1,input_fmt*40
      character cl*1
      integer spec_flag,colindyear,colindday,idum

      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mauxcol ! Avoid compiler warning (unused parameter)
      idum=mcolvav ! Avoid compiler warning (unused parameter)
      idum=mgas    ! Avoid compiler warning (unused parameter)
      idum=mlev    ! Avoid compiler warning (unused parameter)
      idum=mrow_qc ! Avoid compiler warning (unused parameter)
      idum=mspeci  ! Avoid compiler warning (unused parameter)
      idum=mvmode  ! Avoid compiler warning (unused parameter)
      idum=ncell   ! Avoid compiler warning (unused parameter)
      idum=maddln  ! Avoid compiler warning (unused parameter)
      idum=mcharhead! Avoid compiler warning (unused parameter)

      spec_flag=0
      colindyear=0  ! prevent compiler warning (may be used uninitialized)
      colindday =0  ! prevent compiler warning (may be used uninitialized)

      version=
     & ' apply_ghost_correction    Version 1.13    2020-03-12   GCT/JLL'
c     & ' apply_ghost_correction      Version 1.0.1   2012-03-08   GCT'
c     & ' apply_ghost_correction      Version 1.0.0   2011-07-19   NMD'

      call get_ggg_environment(gggdir,dl)

      if (iargc() == 0) then
         write(*,*)'Enter name of input file (e.g. paIn.vav.ada.aia):'
         read(*,'(a)') inputfile
      elseif (iargc() == 1) then
         call getarg(1, inputfile)
      else
         stop 'Usage: $gggpath/bin/apply_ghost_correction aiafile'
      endif
      li=lnbc(inputfile)
      site=inputfile(1:2)
      if(inputfile(li-3:li).ne.'.aia') write(*,*)
     & 'Warning: input file is not of expected type (.aia)'
      outputfile=inputfile(:lnbc(inputfile))//'.gaa'
      open(lunr,file=inputfile, status='old')
      open(lunw,file=outputfile,status='unknown')

c  Open the xx_ghost_corr.dat file and read in the ghost correction
c  factors (GCF) and break times for each gas. 'xx' is the site notation
      filename=gggdir(:lnbc(gggdir))//'tccon/'//site//
     & '_ghost_corr.dat'
      inquire(file=filename,EXIST=fileexists)
      if (.not.fileexists) then
        write(*,*)'The file '//filename(:lnbc(filename))//
     & ' does not exist.'
        stop 
      endif
      open(luns,file=filename,status='old')
      read(luns,*)ncoml,ngas,ngrow
      do k=2,ncoml
         read(luns,*)
      end do
      do k=1,ngrow
         read(luns,*,end=88) endyear(k),endday(k),
     & (gcf(k,kk),gcfe(k,kk),scl(k,kk),kk=1,ngas)
      end do
88    continue
      call substr('xco2 xch4 xn2o xco xh2o',gasname,ngas,ngas1) 
! this is a hard-coded gas order, which should probably be made more flexible/robust
      if (ngas.ne.ngas1) write(*,*)'ngas.ne.ngas1',ngas,ngas1
c  Read the header of the .aia file and figure out the
c  mapping between the gases in the corrections.dat
c  and those in the .vav file header
      read(lunr,countfmt) ncoml,ncol,nrow,naux
      write(lunw,countfmt) ncoml+1+ngrow+1,ncol,nrow,naux
      write(lunw,'(a)') version
      if(ncol.gt.mcol) stop 'increase mcol'
      do j=2,ncoml-2
         read(lunr,'(a)') header
         write(lunw,'(a)') header(:lnbc(header))
      end do
      read(lunr,'(7x,a)') input_fmt
      write(lunw,'(a)') 'format='//input_fmt
      outfmt=input_fmt
c     write(*,*) endyear(1),ngrow
      
      write(lunw,'(a)')'Ghost correction factors:'
      do j=1,ngrow
        write(lunw,'(2i4,5(a,2e11.3))')endyear(j),endday(j),
     & ' xco2',gcf(j,1)*scl(j,1),gcfe(j,1)*scl(j,1),
     & ' xch4',gcf(j,2)*scl(j,2),gcfe(j,2)*scl(j,2),
     & ' xn2o',gcf(j,3)*scl(j,3),gcfe(j,3)*scl(j,3),
     & ' xco', gcf(j,4)*scl(j,4),gcfe(j,4)*scl(j,4),
     & ' xh2o',gcf(j,5)*scl(j,5),gcfe(j,5)*scl(j,5)
      enddo

      read(lunr,'(a)') header
      write(lunw,'(a)') header(:lnbc(header))
      call lowercase(header)
      if (index(header,'spectrum') .gt. 0) spec_flag=1
      call substr(header,headarr,mcol,kcol)
      if(kcol.ne.ncol ) stop 'ncol/kcol mismatch'
      do icol=1,ncol
         do k=1,ngrow
            cf(icol,k)=0.0
         end do
         do kgas=1,ngas
c         write(*,*)kgas,icol,gasname(kgas),headarr(icol),gcf(kgas),
c     &  cf(icol)
           gaserr=gasname(kgas)(:lnbc(gasname(kgas)))//'_error'
           if( headarr(icol) .eq. gasname(kgas) ) then
             do k=1,ngrow
               cf(icol,k)=gcf(k,kgas)*scl(k,kgas)
c               write(*,*) cf(icol,k),gcf(k),k,icol
             enddo
           elseif( headarr(icol) .eq. gaserr ) then
             do k=1,ngrow
               cf(icol,k)=0.0
c               write(*,*) cf(icol,k),gcf(k),k,icol
             enddo
           elseif( headarr(icol) .eq. 'year' ) then
             colindyear=icol
           elseif( headarr(icol) .eq. 'day' ) then
             colindday =icol
           endif
         end do
c        if(icol.gt.naux+spec_flag) write(*,'(i4,1x,a16,f8.3)')
c    &   icol,headarr(icol),cf(icol,k-1)
c         if (k.gt.1) write(*,*) cf(icol,k-1),k
      end do

c  Read each day of data into memory and add the appropriate corrections
c  to the XGas values.
      do irow=1,9999999
         if (spec_flag .eq. 1) then
            read(lunr,*,end=99) specname, (yrow(j),j=1+spec_flag,ncol)
         else
c           write(*,*)'yrow(j),ncol=',yrow,ncol 
            read(lunr,input_fmt,end=99) cl,(yrow(j),j=1,ncol)
         endif
         do k=naux+spec_flag+1,ncol-1,2
           do j=1,ngrow
c             write(*,*) int(yrow(colindyear))
             if (int(yrow(colindyear)).lt.endyear(j)) then
c               write(*,*) endyear(j),cf(k,j),j,k
               yrow(k)=yrow(k)+cf(k,j)
               yrow(k+1)=yrow(k+1)+cf(k+1,j)
c              write(*,*) j,ngrow,cf(k,j)
               goto 97
             elseif (int(yrow(colindyear)).eq.endyear(j)) then
               if (int(yrow(colindday)).le.endday(j)) then
                 yrow(k)=yrow(k)+cf(k,j)
                 yrow(k+1)=yrow(k+1)+cf(k+1,j)
c                 write(*,*)'b',j,ngrow,cf(k,j)
                 goto 97
               endif
             endif
           enddo
97         continue !write(*,*) j,ngrow,cf(k,j)
         end do
         if (spec_flag .eq. 1) then
            write(lunw,outfmt) specname,(yrow(j),j=1+spec_flag,ncol)
         else
            write(lunw,outfmt) cl,(yrow(j),j=1,ncol)
         endif
      end do         ! do irow=1,9999999
c      write(*,*) cf
      stop ' irow exceeded 9999999'
99    close (lunr)
      close(lunw)
      stop
      end
