C-------------------------------------------------
C  Merges the Hartmann/Tran/Niro codes for CO2, O2 and CH4 Line Mixing.
c
c  Program asks for the .ggg file for the window under consideration.
c  From this it computes the range of frequencies over which the
c  absorption coefficients must be pre-computed.
c
c  It also reads the .mav file to get the P/T information for each level.
c  
C-------------------------------------------------
C
C     INPUT VARIABLES 
C --------------
C sigmin [cm-1]
C sigmax [cm-1]
C xCO2 [no unit]
C T [K]
C p [atm]
C  
C     MixFull = Switch to full diagonalization line-mixing
C     rdmmult = Distance from line center at which you
C               can shift Voigt to lorentz (save CPU time)
C               values lower than 30 should be avoided in
C               order to minimize error
C
C   RESULTS 
c     AbsV:  Absorption Coefficient assuming Voigt lineshape
c            (neglecting LM)
c
c     AbsY:  Difference between 1st order LM and Voigt calculations
c
c     AbsW: CO2: Difference between Full and 1st order LM
c           CH4: Voigt Calculation using Frankenberg linelist
c            O2: CIA absorption coefficients
c
c  The results (AbsV, AbsY, AbsW) are written to a binary output file
C
      implicit none
      include "../gfit/ggg_int_params.f"
      integer lunr,lunr_mav,lun_vac,lun_wbin,
     & nlhead,mcol,ncol,kcol,icol,nn,
     & nlev_mav,ilev,k,nSigmx,kcol_temp,kcol_pres,
     & lcolon,k1,k2,fbc,fnbc,lnbc,lrt,irec,
     & kcol_h2o, kcol_co2, kcol_ch4,kcol_o2
      parameter(lunr=13,lunr_mav=14,lun_vac=15,lun_wbin=16,
     & nSigmx=500000,mcol=200)
      logical MixFull
      real*8 StotMax,vv(mcol)
C Results (Absorption Coefficients)
      real*8 AbsV(nsigmx)
      real*8 AbsY(nsigmx)
      real*8 AbsW(nsigmx)
      real*8 sigmin, sigmax, dsig, sig, resmax, frqcen, width
      parameter (resmax=0.625d0)
      integer i,j,lmav,nsig,kcp1,kcp2,nsh,nhwmax,nscycle,idum
      parameter (nscycle=25)
      character mavfile*80, outfile*80, next_spectrum*12,
     & version*80,
     & gggdir*(mpath),dl*1,
     & path_to_input_files*(mpath+80),
     & input_ggg*80, winfo*200, gas*12,
     & header_string*2200, header_vector(mcol)*12,
     & path_to_bandinfo*(mpath+80),
     & path_to_wfile*(mpath+80),
     & path_to_sfile*(mpath+80)

      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mauxcol ! Avoid compiler warning (unused parameter)
      idum=mcolvav ! Avoid compiler warning (unused parameter)
      idum=mgas    ! Avoid compiler warning (unused parameter)
      idum=mlev    ! Avoid compiler warning (unused parameter)
      idum=mrow_qc ! Avoid compiler warning (unused parameter)
      idum=mspeci  ! Avoid compiler warning (unused parameter)
      idum=mvmode  ! Avoid compiler warning (unused parameter)
      idum=ncell   ! Avoid compiler warning (unused parameter)
      idum=nchar   ! Avoid compiler warning (unused parameter)

      version=' compute_lm_absco    Version 1.11    2018-09-12   GCT'
      write(*,*) version
      MixFull=.false.
      MixFull=.true.

c     Platform specification:
      call get_ggg_environment(gggdir, dl)
      lrt=lnbc(gggdir)     !Length of gggdir
c
C----------
c
c T and P vertical profile
c
      if (iargc() == 0) then
         write(*,'(a)')' Enter name of input .ggg file:'
         read(*,'(a)') input_ggg
      elseif (iargc() == 1) then
         call getarg(1, input_ggg)
      else
         stop 'Usage: $gggpath/bin/compute_lm_absco gggfile'
      endif
      open(lunr,file=input_ggg, status='old')
      read(lunr,*)nlhead
      do j=2,nlhead
        read(lunr,'(a)') winfo
        lmav=index(winfo,'.mav')
        if(lmav.gt.0) mavfile=winfo(:lmav+4)
      end do
      write(*,*) winfo
      if(index(winfo,' lm ').eq.0)then
        stop ' The "lm" switch was not found in the ggg file.'
      else
        write(*,*)' Found the "lm" switch in the ggg file. Proceeding.'
      endif

      nn=0
      lcolon=index(winfo,':')
88    k1=fnbc(winfo(lcolon+1:))+lcolon
      k2=fbc(winfo(k1+1:))+k1-1
      if(winfo(k1:k2).eq.'1co2') k1=k1+1  ! drop the prefix 1
      if(winfo(k1:k2).eq.'1ch4') k1=k1+1  ! drop the prefix 1
      if(winfo(k1:k2).eq.'1o2')  k1=k1+1  ! drop the prefix 1
      if(winfo(k1:k2).eq.'ao2')  k1=k1+1  ! drop the prefix a
      if(winfo(k1:k2).eq.'bo2')  k1=k1+1  ! drop the prefix b
      if(winfo(k1:k2).eq.'wco2') k1=k1+1   ! drop the prefix w
      if(winfo(k1:k2).eq.'lco2') k1=k1+1   ! drop the prefix l
      gas=winfo(k1:k2)
      if(gas.ne.'co2' .and. gas.ne.'ch4' .and. gas.ne.'o2') then
         lcolon=k2
         nn=nn+1
         if(nn.gt.5) stop 'Unrecognised gas: nn>5'
         go to 88
      endif
c      write(*,*)'gas=',gas
      
c      write(*,*) winfo(k1:k2)
c      write(*,*)lcolon,k1,k2
      path_to_input_files=gggdir(:lrt)//'src/compute_lm_absco/'
     &//'data_'//winfo(k1:k2)//'/'
c      path_to_input_files='/home/toon/hartmann/'//
c     & winfo(k1:k2)//'/data_'//winfo(k1:k2)//'/'

      path_to_bandinfo=path_to_input_files
      path_to_wfile=path_to_input_files
      path_to_sfile=path_to_input_files

c      write(*,*)path_to_input_files
      read(winfo,*) frqcen,width
      
      dsig=0.666666d-06*frqcen
      kcp1=int((frqcen-width/2)/dsig)
      kcp2=int((frqcen+width/2)/dsig)
      nsh=int(2+2*resmax/dsig)
      nhwmax=nint(nscycle*resmax/dsig)
      nsig=kcp2-kcp1+2*nhwmax+2*nsh
      write(*,*)'nus,nue,grid,ncp=',frqcen-width/2,frqcen+width/2,
     &  dsig,nsig,nhwmax,nsh


      sigmin=dsig*(kcp1-nsh-nhwmax+1)
      sigmax=dsig*(kcp1-nsh-nhwmax+nsig)
c      write(*,*) sigmin, sigmax, nsig, dsig

c--------- 
      StotMax=0.4E-21  ! The weak CO2 bands at 6220 & 6338 cm-1 are 0.44E-21
      if(gas.eq.'co2') then
         call DetBand(path_to_bandinfo,sigmin,sigmax,StotMax)
         call ReadW(path_to_wfile)
      endif
c
      open(lunr_mav,file=mavfile,status='old')
      read(lunr_mav,*)
1     read(lunr_mav,'(14x,a12)',end=99) next_spectrum
      read(lunr_mav,*) nlhead,ncol,nlev_mav
      do k=2,nlhead
         read(lunr_mav,'(a)') header_string
      end do
      call substr(header_string, header_vector, mcol,kcol)
      if(kcol.ne.ncol) then
         write(*,*) 'kcol, ncol =', kcol,ncol
           stop 'kcol/ncol mismatch'
      endif
      kcol_h2o=0
      kcol_co2=0
      kcol_ch4=0
      kcol_o2=0
      kcol_pres=0
      kcol_temp=0
      do icol=1,ncol
        if(index(header_vector(icol),'Temp ').gt.0) kcol_temp=icol
        if(index(header_vector(icol),'temp ').gt.0) kcol_temp=icol
        if(index(header_vector(icol),'Pres ').gt.0) kcol_pres=icol
        if(index(header_vector(icol),'pres ').gt.0) kcol_pres=icol
        if(index(header_vector(icol),'1h2o ').gt.0) kcol_h2o=icol
        if(index(header_vector(icol),'1co2 ').gt.0) kcol_co2=icol
        if(index(header_vector(icol),'1ch4 ').gt.0) kcol_ch4=icol
        if(index(header_vector(icol), '1o2 ').gt.0) kcol_o2 =icol
      end do

      write(outfile,'(a12,a1,a,a1,i5.5,a1,i5.5,a3)')  next_spectrum,'_',
     & winfo(k1:k2),'_',nint(sigmin),'_',nint(sigmax),'_lm'
      open (unit=lun_vac,file='compute_lm_absco.rpt',
     & status='unknown')
      write(lun_vac,*)2,7
      write(lun_vac,*)' i t p f v y w'
c  Open binary absco file
      open(lun_wbin,file=outfile(:lnbc(outfile))//'.bin',
     &access='direct',status='unknown',form='unformatted',recl=12*nsig)
      irec=0
      do ilev=1,nlev_mav
         read(lunr_mav,*) (vv(icol),icol=1,ncol)
         if(vv(kcol_pres).le.0.0) vv(kcol_pres)=0.0001
         write(*,'(i4, 2f9.1, 2f11.8)') ilev, vv(1), vv(2), vv(3)
         if(gas.eq.'o2'.or.gas.eq.'ao2') then
            call LMandCIAO2(path_to_input_files,
     &      vv(kcol_temp),vv(kcol_pres),vv(kcol_o2),vv(kcol_h2o),
     &      sigmin,dSig,nSig,AbsV, AbsY, AbsW)
         elseif(gas.eq.'co2' .or. gas.eq.'1co2') then
            call CompAbsCO2(path_to_sfile,
     &      vv(kcol_temp),vv(kcol_pres),vv(kcol_co2),vv(kcol_h2o),
     &      sigmin,dsig,nsig,AbsV,AbsY,AbsW,MixFull)
c   AbsV is the Voigt Lineshape absorption coeffs
c   AbsY is the difference: 1'st order LM - Voigt abscoeffs
c   AbsW is the difference: Full LM  -  1'st order LM abscoeffs (CO2)
c   AbsW is the CIA: (O2)
c   AbsW is the non-LM absorption (CH4)
c   So  AbsY+AbsW = Full LM - Voigt abscoeffs
         elseif(gas.eq.'ch4') then
c            write(*,*) 'calling CompAbsCH4...'
            call CompAbsCH4(path_to_input_files,sigmin,dsig,nsig,
     &      vv(kcol_temp),vv(kcol_pres),vv(kcol_ch4),vv(kcol_h2o),
     &      AbsV, AbsY, AbsW)
         else
            write(*,*)'No LM subroutine for: ',gas
            stop
         endif

         write(lun_wbin,rec=ilev)
     &   (sngl(AbsV(i)),sngl(AbsY(i)),sngl(AbsW(i)),i=1,nsig)
         do i=1,nsig
            sig=sigmin+(i-1)*dsig
         write(lun_vac,123)ilev,vv(2),vv(3),sig,AbsV(i),AbsY(i),AbsW(i)
 123     format(i4,f8.2,f8.5,f12.5,3(1pe12.3))
         end do  ! i=1,nsig
         irec=irec+nsig
      end do   ! ilev=1,nlev_mav
      close(lun_vac)
      go to 1   !  Loop to next "Next Spectrum"
99    close(lunr_mav)
      close(lun_wbin)
      stop
      end
