c  Program to create a GGG-compatible sunrun from two inputs:
c  1) An ascii list of OPUS spectra 
c  2) A site-dependent file containing the invarient parameter (e.g. lat, long, alt,..)
c
c  To create the ascii input list, use the program list_maker.
c
c  The sunrun program can then be immediately run from there.
c  It will search the local directory for the .gnd file, and
c  write the .gop file to the sunruns/ directory.
c
      implicit none
      include "../ggg_int_params.f"
      include "params.f"

      integer*4 iy,im,id,hh,mm,ss,ms,pkl,prl,gfw,gbw,
     & fnbc,lnbc,fbc,ispe,iend,dtype,nsp,nip,dfr,bytepw,
     & mdet,ndet,idet,             ! Number of detectors
     & ifirst,ilast,
     & lr,lrt,ls,lunr,doy,
     & possp,istat,object,mcol,ncol,instr,lst
      parameter (lunr=14,mcol=40,mdet=4)
c
      real*8 tins,pins,hins,tout,pout,pout_corr,hout,
     & wspd,wdir,gmt,fovi,opd,snr,asza,delwav,
     & wavtkr,oblat,oblon,obalt,lfl,hfl,foc,
     & fsf,tcorr,vdc,lse,lsu,
     & nus(mdet),nue(mdet),
     & lasf,sia,sis,fvsi,aipl,tel_mag,
     & fxv,lxv,apt,dur,vel,phr,res,ptrue
c
      integer k1(mdet), k2(mdet),irec, nhead, lsss

      character 
     & sss(mdet)*2,
     & parfile*(mfilepath),
     & header*512,outarr(mcol)*20
      character
     & col1*1,                    !first column of runlog record
     & apf*2,                     !apodization function (e.g. BX N2, etc)
     & dl*1,                      !forward or backward slash
     & ext*3,                     !geometry ['air','bal','gnd','lab',orb','syn']
     & gggdir*(mpath),            !ggg directory path (GGGPATH?)
     & specname*(nchar),          !spectrum name
     & version*64                 !current program version

c
      version=' create_sunrun     Version 2.1.2     12-Apr-2012     GCT'
      write(6,'(a)') version

      col1=' '
      call getendian(iend)  ! iend=+1 on Sun; iend=-1 on PC
      vdc =0.d0  ! initialize VDC

c  Prompt user for name of input list
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

      call get_ggg_environment(gggdir,dl)
      lrt=lnbc(gggdir)       !Length of root
c
c  Read file containing the invariant parameters.
      parfile=gggdir(:lrt)//'tccon/'//logfile(1:2)//'_sunrun.dat'
      write(*,*) parfile
      open(luns,file=parfile,status='old')
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cc     Root specification:        DG000909
c      call getenv('GGGPATH',root)
c      dl='/'
c      root=root(:lnbc(root))//dl
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      header=
     &' Spectrum_File_Name                     Obj   tcorr'//
     &'   oblat   oblon   obalt   tins   pins   hins'//
     &'  tout  pout   hout    sia    fvsi   wspd   wdir'//
     &'   Nus    Nue      FSF      lasf     wavtkr'//
     &'   AIPL   TM'
      call substr(header, outarr, mcol,ncol)

      read(luns,*) nhead
      do irec=2,nhead
         read(luns,*)
      end do
      read(luns,*) instr
      read(luns,*) object
      read(luns,*) tcorr
      read(luns,*) oblat
      read(luns,*) oblon
      read(luns,*) obalt
      read(luns,*) tins
      read(luns,*) pins
      read(luns,*) hins
      read(luns,*) tout
      read(luns,*) pout
      read(luns,*) pout_corr
      read(luns,*) hout
      read(luns,*) sia
      read(luns,*) fvsi
      read(luns,*) wspd
      read(luns,*) wdir
      read(luns,*) ndet
      do idet=1,ndet
         read(luns,*) k1(idet),k2(idet),sss(idet),nus(idet),nue(idet)
      end do
      read(luns,*) fsf
      read(luns,*) lasf
      read(luns,*) wavtkr
      read(luns,*) aipl
      read(luns,*) tel_mag 
      close(luns)
      write(*,*) 'Closed xx_sunrun.dat file'
      lsss=lnbc(sss(1))

c      outfile=gggdir(:lrt)//'sunruns'//dl//ext//dl//logfile(:lr-2)//'m4'
      outfile=gggdir(:lrt)//'sunruns'//dl//ext//dl//logfile(:lr-2)//'op'
      open(lunt,file=outfile,status='unknown')
      write(lunt,*)3,ncol
      write(lunt,'(a)') version
      write(lunt,'(a)')header(:lnbc(header))
c
      open(lunr,file=logfile,status='old')
      do ispe=1,999999  !---------Main loop over spectra----------
         read(lunr,'(a)',end=99) string 
c         write(*,*) string(:lnbc(string))

c  Remove any leading spaces/blanks.
         string=string(fnbc(string):)
         ls=fbc(string)
         specname=string(:ls-1)
         do idet=1,ndet
            if(specname(k1(idet):k2(idet)).eq.sss(idet)(1:lsss)) exit
         end do
c
c  find the spectral file, return the PATH to the spectrum
          dplist=gggdir(:lrt)//'config'//dl//'data_part.lst'
          call gindfile(dplist,specname,path)
          if(lnbc(path).le.0) then
             write(6,*) ' Not Found : ',specname
             stop
          endif

c          if(logfile(1:2).eq.'m4') then
          if(instr.eq.1) then
            call read_mkiv_header(specname,path,iend,ifirst,ilast,possp,
     &      bytepw,apf,delwav,opd,fovi,snr,oblat,oblon,obalt,pout,
     &      tout,hout,asza,iy,im,id,gmt,wavtkr,tins,pins,hins,lasf)
          else
            call read_opus_header(path,iend,dtype,nsp,fxv,lxv,iy,im,
     &      id,hh,mm,ss,ms,apt,dur,vel,apf,phr,res,lasf,foc,nip,dfr,
     &      pkl,prl,gfw,gbw,lfl,hfl,possp,oblat,oblon,obalt,
     &      tins,pins,hins,tout,pout,hout,wspd,wdir,sia,sis,vdc,
     &      lst,lse,lsu)
          endif

c  Apply correction to measured surface pressure.
          ptrue = pout + pout_corr

c  Calculate fvsi as sis/sia if sia is not zero or missing, and if
c  it's not a lamp run:
          if(sia.le.0.0 .or. ext(:1).eq.'l') then
            if(vdc.gt.0) then
              fvsi=vdc
            else
c             Use the missing value in xx_sunrun.dat
            endif
          else
            fvsi=sis/sia
          endif
          if(fvsi.ge.1.0)fvsi=1.
c  moving it here assumes the default value is set to something <1

c  If the TCCON spectrum has no red filter, limit the upper
c  frequency (nue) to 13500.0d0
          if(instr.eq.2) then ! TCCON spectrum
             if(specname(15:15).eq.'0') then ! no Si red filter
                  nue=min(nue,13500.0d0)
             endif
          endif

c  Site-dependent corrections should go in here. For example, to correct 
c  a timing error in Lamont:

c  Apply correction for timing errors from July 15 - October 24, 2010
          if(instr.eq.2 .and. logfile(1:2).eq.'oc') then ! Lamont
            if(iy.eq.2010 .and. ((im.eq.7 .and. id.ge.15) .or.
     &        im.eq.8 .or. im.eq.9 .or.
     &        (im.eq.10 .and. id.le.24))) then
c              write(*,*)'Applying Lamont timing error correction.'
c              doy 181 is June 30, 196 is July 15, 297 is October 24
c              31 is approximate length of month in days
               doy = 181 + (im-7)*31 + id
c              linearly decreasing tcorr from 0 to -37 seconds
               tcorr = (doy-196)*37/(196-297)
            else
               tcorr = 0.0d0
            endif
          endif


          call write_sunrun(lunt,col1,specname,object,tcorr,oblat,
     &    oblon,obalt,tins,pins,hins,tout,ptrue,hout,sia,fvsi,wspd,wdir,
     &    nus(idet),nue(idet),fsf,lasf,wavtkr,aipl,tel_mag,istat)
c
      if(mod(ispe,1000).eq.0) write(*,*)ispe
      end do ! -------------Main loop over spectra----------------
c
 99   close(lunr)
      close(lunt)
      write(*,*) ispe-1, ' spectra found'
      stop
      end
