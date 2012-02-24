c  Program to create a GGG-compatible sunrun from an ascii list
c  of spectra.
c
      implicit none
      include "../../ggg_int_params.f"
      include "../params.f"

      integer*4 iy,im,id,ifirst,ilast,bytepw,fnbc,lnbc,fbc,ispe,iend,
     & lr,lrt,ls,
     & possp,istat,object
c
      real*8 tins,pins,hins,tout,pout,hout,
     & oblat,oblon,obalt,wavtkr,gmt,delwav,opd,fovi,snr,asza,
     & fsf,tcorr,nus,nue,lasf,sia,sis,aipl,tel_mag
c
      character 
     & col1*1,                    !first column of runlog record
     & apf*2,                     !apodization function (e.g. BX N2, etc)
     & dl*1,                      !forward or backward slash
     & ext*3,                     !geometry ['air','bal','gnd','lab',orb','syn']
     & gggdir*(mpath),            !ggg directory path (GGGPATH?)
     & specname*(nchar)           !spectrum name

c
      write(6,*)'SUNRUN     Version 1.1.0    31-Oct-2005    GCT'
      col1=' '
c
      call get_ggg_environment(gggdir, dl)
      lrt=lnbc(gggdir)       !Length of root
      call getendian(iend)
      lr=0
      do while(lr.eq.0)
         write(6,'(a)') 'Enter name of input file (e.g. pa2004.gnd): '
         read(*,'(a)') logfile
         lr=lnbc(logfile)
      end do
      ext=logfile(lr-2:lr)

      open(luns,file=logfile,status='old')
      outfile=gggdir(:lrt)//'sunruns'//dl//ext//dl//logfile(:lr-2)//'m4'
      open(lunt,file=outfile,status='unknown')
      write(lunt,'(a)') ' Spectrum_File_Name   Obj   tcorr    oblat'//
     &'   oblon   obalt   tins  pins   hins   tout  pout   hout'//
     &'    Nus    Nue      FSF       lasf     wavtkr   sia   sis'//
     &'  aipl   tm'
c
      object=2 ! Sun
      tcorr=0.0d0       ! Times are in UT
      fsf=1.0d0
      sia=0.0d0
      sis=0.0d0
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

          call read_mkiv_header (specname,path,iend,ifirst,ilast,possp,
     &    bytepw,apf,delwav,opd,fovi,snr,oblat,oblon,obalt,
     &    pout,tout,hout,asza,iy,im,id,gmt,wavtkr,tins,pins,hins,lasf)

      aipl=0.001        ! Airmass-Independent Path Length (km)

          call write_sunrun(lunt,col1,specname,object,tcorr,oblat,
     &    oblon,obalt,tins,pins,hins,tout,pout,hout,nus,nue,fsf,
     &    lasf,wavtkr,sia,sis,aipl,tel_mag,istat)
c
      end do ! -------------Main loop over spectra----------------
c
 99   close(luns)
      close(lunt)
      write(*,*) ispe-1, ' spectra found'
      stop
      end
