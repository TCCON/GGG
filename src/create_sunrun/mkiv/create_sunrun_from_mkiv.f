c  Program to create a GGG-compatible sunrun from an ascii list
c  of spectra.
c
      implicit none
      include "../../gfit/ggg_int_params.f"
      include "../params.f"

      integer*4 iy,im,id,ifirst,ilast,bytepw,fnbc,lnbc,fbc,ispe,iend,
     & lr,lrt,ls,mcol,ncol,
     & possp,istat,object,idum
      parameter (mcol=40)
c
      real*8 tins,pins,hins,tout,pout,hout,wspd,wdir,
     & oblat,oblon,obalt,wavtkr,gmt,delwav,opd,fovi,snr,asza,
     & fsf,tcorr,nus,nue,lasf,sia,fvsi,aipl,tel_mag
c
      character 
     & col1*1,                    !first column of runlog record
     & apf*2,                     !apodization function (e.g. BX N2, etc)
     & dl*1,                      !forward or backward slash
     & ext*3,                     !geometry ['air','bal','gnd','lab',orb','syn']
     & version*60,
     & header*512,outarr(mcol)*20,
     & gggdir*(mpath),            !ggg directory path (GGGPATH?)
     & specname*(nchar)           !spectrum name

      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mauxcol  ! Avoid compiler warning (unused parameter)
      idum=mcolvav  ! Avoid compiler warning (unused parameter)
      idum=mgas     ! Avoid compiler warning (unused parameter)
      idum=mlev     ! Avoid compiler warning (unused parameter)
      idum=mrow_qc  ! Avoid compiler warning (unused parameter)
      idum=mspeci   ! Avoid compiler warning (unused parameter)
      idum=mvmode   ! Avoid compiler warning (unused parameter)
      idum=ncell    ! Avoid compiler warning (unused parameter)
      idum=nchar    ! Avoid compiler warning (unused parameter)

      version=
     & 'create_sunrun_from_mkiv   Version 1.2.0    01-Aug-2012   GCT'
      col1=' '
c
      call get_ggg_environment(gggdir, dl)
      lrt=lnbc(gggdir)       !Length of root
      call getendian(iend)
      lr=0
      do while(lr.eq.0)
         if (iargc() == 0) then
            write(6,'(a)') 'Enter name of input file (e.g. pa2004.gnd):'
            read(*,'(a)') logfile
         elseif (iargc() == 1) then
            call getarg(1, logfile)
         else
            stop 'Usage: $gggpath/bin/create_sunrun pa2004.gnd'
         endif

         lr=lnbc(logfile)
      end do
      ext=logfile(lr-2:lr)

      header= ' Spectrum_File_Name                        '//
     &'                Obj   tcorr    oblat  oblon   obalt'//
     &'   tins  pins   hins'//
     &'   tout  pout   hout   sia   fvsi   wspd   wdir    Nus    Nue'//
     &'      FSF       lasf     wavtkr   AIPL   TM'
      call substr(header, outarr, mcol,ncol)

      open(luns,file=logfile,status='old')
      outfile=gggdir(:lrt)//'sunruns'//dl//ext//dl//logfile(:lr-2)//'m4'
      open(lunt,file=outfile,status='unknown')
      write(lunt,*)3,ncol
      write(lunt,'(a)') version
      write(lunt,'(a)')header(:lnbc(header))

c
      object=2 ! Sun
      tcorr=0.0d0       ! Times are in UT
      fsf=1.0d0
      sia=0.0d0
      fvsi=0.0d0
      tel_mag=1.0d0
      wavtkr=9900.0d0
c
      do ispe=1,999999  !---------Main loop over spectra----------
         read(luns,'(a)',end=99) string 

c  Remove any leading spaces/blanks.
         string=string(fnbc(string):)
         ls=fbc(string)
         specname=string(:ls-1)
         if(specname(2:3).eq.'hg') then      ! HgCdTe detector
            nus=650.0d0
            nue=1850.0d0
         elseif(specname(2:3).eq.'in') then  ! InSb detector
            nus=1850.0d0
            nue=5650.0d0
         else
            write(*,*) 'Unknown spectrum type: '//string
            stop
         endif
c
c  find the spectral file, return the PATH to the spectrum
          dplist=gggdir(:lrt)//'config'//dl//'data_part.lst'
          call gindfile(dplist,specname,path)
          if(lnbc(path).le.0) then
            write(6,*) ' Not Found : ',specname
            stop
          endif

           write(*,*)'Read_mkiv header'
          call read_mkiv_header (specname,path,iend,ifirst,ilast,possp,
     &    bytepw,apf,delwav,opd,fovi,snr,oblat,oblon,obalt,
     &    pout,tout,hout,asza,iy,im,id,gmt,wavtkr,tins,pins,hins,lasf)
          write(*,*) specname, iy, im, id, gmt,asza

      aipl=0.001        ! Airmass-Independent Path Length (km)

           write(*,*)'write_sunrun'
          call write_sunrun(lunt,col1,specname,object,tcorr,oblat,
     &    oblon,obalt,tins,pins,hins,tout,pout,hout,sia,fvsi,
     &    wspd,wdir,nus,nue,fsf,lasf,wavtkr,aipl,tel_mag,istat)


      end do ! -------------Main loop over spectra----------------
c
 99   close(luns)
      close(lunt)
      write(*,*) ispe-1, ' spectra found'
      stop
      end
