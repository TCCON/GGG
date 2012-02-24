c  Program to create a GGG-compatible sunrun from an ascii list
c  of OPUS spectra. Contains corrections specific to Parkfalls.
c
c  To create the input file for this program run list_maker.
c
c  The sunrun program can then be immediately run from there.
c  It will search the local directory for the .gnd file, and
c  write the .gop file to the sunruns/ directory.
c
      implicit none
      include "../../ggg_int_params.f"
      include "../params.f"

      integer*4 iy,im,id,hh,mm,ss,ms,pkl,prl,gfw,gbw,
     & fnbc,lnbc,fbc,ispe,iend,dtype,nsp,nip,dfr,
     & lr,lrt,ls,
     & possp,istat,object
c
      real*8 tins,pins,hins,tout,pout,hout,
     & wavtkr,oblat,oblon,obalt,lfl,hfl,foc,
     & fsf,tcorr,nus,nue,lwn,sia,sis,aipl,tel_mag,
     & fxv,lxv,apt,dur,vel,phr,res,pout_cor
c
      character 
     & spfmt*2,
c    & user*8,
     & col1*1,                    !first column of runlog record
     & apf*2,                     !apodization function (e.g. BX N2, etc)
     & dl*1,                      !forward or backward slash
     & ext*3,                     !geometry ['air','bal','gnd','lab',orb','syn']
     & gggdir*(mpath),            !ggg directory path (GGGPATH?)
     & specname*(nchar)           !spectrum name

c
       pout_cor = 0.0

      write(6,*)'SUNRUN     Version 1.1.0    31-Oct-2005    GCT'
      write(6,'(a,f5.2,a)')'Pout will be adjusted by ',pout_cor,' mb'
      col1=' '
      call getendian(iend)  ! iend=+1 on Sun; iend=-1 on PC
c
      call get_ggg_environment(gggdir, dl)
      lrt=lnbc(gggdir)       !Length of root
      lr=0
      do while(lr.eq.0)
         write(6,'(a)') 'Enter name of input file (e.g. pa2004.gnd): '
         read(*,'(a)') logfile
         lr=lnbc(logfile)
      end do
      ext=logfile(lr-2:lr)
      spfmt=logfile(lr-1:lr)

      open(luns,file=logfile,status='old')
      outfile=gggdir(:lrt)//'sunruns'//dl//ext//dl//logfile(:lr-2)//'op'
      open(lunt,file=outfile,status='unknown')
      write(lunt,'(a)') ' Spectrum_File_Name   Obj   tcorr    oblat'//
     &'   oblon   obalt   tins  pins   hins   tout  pout   hout'//
     &'    Nus    Nue      FSF       lasf     wavtkr   sia   sis'//
     &'   AIPL   TM'
c
      object=2 ! Sun
      tcorr=0.0d0       ! Times are in UT
      fsf=0.99999999d0
      wavtkr=9900.0d0
      aipl=0.002        ! Airmass-Independent Path Length (km)
      tel_mag=1.0
c
      do ispe=1,999999  !---------Main loop over spectra----------
         read(luns,'(a)',end=99) string 

c  Remove any leading spaces/blanks.
         string=string(fnbc(string):)
         ls=fbc(string)
         specname=string(:ls-1)
         if(specname(16:16).eq.'a') then      ! InGaAs detector
            nus=4000.0d0
            nue=11000.0d0
         elseif(specname(16:16).eq.'b') then  ! Si detector
            if(specname(15:15).eq.'0') then
               nus=11000.0d0  ! No Red Filter
               nue=13500.0d0  ! No Red Filter
            else
               nus=11000.0d0  ! Red Filter
               nue=15000.0d0  ! Red Filter
            endif
         elseif(specname(16:16).eq.'c') then   ! InSb detector
            nus=1800.0d0
            nue=6000.0d0
         else
            write(*,*) 'Spectrum Name Length =',ls-1
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

       call read_opus_header(path,iend,dtype,nsp,fxv,lxv,iy,im,
     & id,hh,mm,ss,ms,apt,dur,vel,apf,phr,res,lwn,foc,nip,dfr,
     & pkl,prl,gfw,gbw,lfl,hfl,possp,oblat,oblon,obalt,
     & tins,pins,hins,tout,pout,hout,sia,sis)

c  Apply correction to measured surface pressure.
          pout = pout + pout_cor

          call write_sunrun(lunt,col1,specname,object,tcorr,oblat,
     &    oblon,obalt,tins,pins,hins,tout,pout,hout,nus,nue,fsf,
     &    lwn,wavtkr,sia,sis,aipl,tel_mag,istat)
c
      if(mod(ispe,1000).eq.0) write(*,*)ispe
      end do ! -------------Main loop over spectra----------------
c
 99   close(luns)
      close(lunt)
      write(*,*) ispe-1, ' spectra found'
      stop
      end
