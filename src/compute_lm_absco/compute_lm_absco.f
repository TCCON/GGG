C-------------------------------------------------
C  Merges the Hartmann/Tran codes for CO2, O2 and CH4 Line Mixing into a single program.
c
c  Program asks for the .ggg file for the window under consideration.
c  From this it computes the range of frequencies over which the
c  absorption coefficients must be pre-computed.
c
c  It also reads the .mav file to get the P/T information for each level.
c  
C-------------------------------------------------
C
C	INPUT VARIABLES 
C --------------
C	 sigmin [cm-1]
C	 sigmax [cm-1]
C	 xCO2 [no unit]
C	 T [K]
C	 p [atm]
C  
C     MixFull = Switch to full diagonalization line-mixing
C     rdmmult = Distance from line center at which you
C               can shift Voigt to lorentz (save CPU time)
C               values lower than 30 should be avoided in
C               order to minimize error
C
C	RESULTS 
C --------------
C     AbsV : Absorption Coefficient neglecting LineMixing
C            (assuming Voigt Line-Shapes) (Cm-1)
C     AbsY : Absorption Coefficient predicted using the First
C            Order Line-Mixing Approximation (Cm-1)
C     AbsW : Absorption Coefficient predicted using Full
C            diagonalization Line-Mixing (Cm-1)
C-------------------------------------------------
c     AbsV:  Contains Voigt lineshape
c
c     AbsY:  Difference between 1st order LM and Voigt calculationa
c
c     AbsW: CO2: Difference between Full and 1st order LM
c           CH4: Voigt Calculation using Frankenberg linelist
c           aO2: CIA absorption coefficients
C
      implicit none
      include "../ggg_int_params.f"
      integer lunr,lun_mav,lun_vac,lun_wbin,
     & nlhead,mcol,ncol,kcol,icol,
     & nlev,ilev,k,Li,nSigmx,kcol_temp,kcol_pres,
     & lcolon,k1,k2,fbc,fnbc,lnbc,lrt,irec,
     & kcol_h2o, kcol_co2, kcol_ch4,kcol_o2
      parameter(lunr=13,lun_mav=14,lun_vac=15,lun_wbin=16,
     & nSigmx=500000,mcol=200)
      logical MixFull
      real*8 pshift,StotMax,vv(mcol)
      real*8 z,temp,pres,dens
C Results (Absorption Coefficients)
      real*8 AbsV(nsigmx)
      real*8 AbsY(nsigmx)
      real*8 AbsW(nsigmx)
      real*8 sigmin, sigmax, dsig, sig, resmax, frqcen, width
      parameter (resmax=0.375d0)
      integer i,j,lmav,nsig,kcp1,kcp2,nsh,nhwmax,nscycle
      parameter (nscycle=25)
      character mavfile*80, outfile*80, next_spectrum*10,
     & version*80,
     & gggdir*(mpath),dl*1,
     & path_to_input_files*(mpath+80),
     & input_ggg*80, winfo*200, gas*12,
     & header_string*2000, header_vector(mcol)*12,
     & path_to_bandinfo*(mpath+80),
     & path_to_wfile*(mpath+80),
     & path_to_sfile*(mpath+80)

      version=' compute_lm_absco     V1.0.3    GCT  2012-02-22'
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
      write(*,'(a)')' Enter name of input .ggg file:'
      read(*,'(a)') input_ggg
      open(lunr,file=input_ggg, status='old')
      read(lunr,*)nlhead
      do j=2,nlhead
        read(lunr,'(a)') winfo
        lmav=index(winfo,'.mav')
        if(lmav.gt.0) mavfile=winfo(:lmav+4)
      end do
      write(*,*) winfo
      lcolon=index(winfo,':')
      k1=fnbc(winfo(lcolon+1:))+lcolon
      k2=fbc(winfo(k1+1:))+k1-1
      if(winfo(k1:k2).eq.'1co2') k1=k1+1  ! drop the prefix 1
      gas=winfo(k1:k2)
      
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
      open(lun_mav,file=mavfile,status='old')
      read(lun_mav,*)
1     read(lun_mav,'(14x,a10)',end=99) next_spectrum
      read(lun_mav,*) nlhead,ncol,nlev
      do k=2,nlhead
         read(lun_mav,'(a)') header_string
      end do
      call substr(header_string, header_vector, mcol,kcol)
      if(kcol.ne.ncol) then
         write(*,*) 'kcol, ncol =', kcol,ncol
           stop 'kcol/ncol mismatch'
      endif
      kcol_co2=0
      kcol_ch4=0
      kcol_o2=0
      do icol=1,ncol
        if(index(header_vector(icol),'Temp').gt.0) kcol_temp=icol
        if(index(header_vector(icol),'Pres').gt.0) kcol_pres=icol
        if(index(header_vector(icol),'1h2o').gt.0) kcol_h2o=icol
        if(index(header_vector(icol),'1co2').gt.0) kcol_co2=icol
        if(index(header_vector(icol),'1ch4').gt.0) kcol_ch4=icol
        if(index(header_vector(icol), '1o2').gt.0) kcol_o2 =icol
      end do

      write(outfile,'(a10,a1,a,a1,i5.5,a1,i5.5,a3)')  next_spectrum,'_',
     & winfo(k1:k2),'_',nint(sigmin),'_',nint(sigmax),'_lm'
      open (unit=lun_vac,file='compute_lm_absco.rpt',
     & status='unknown')
      write(lun_vac,*)2,7
      write(lun_vac,*)' i t p f v y w'
c  Open binary absco file
      open(lun_wbin,file=outfile(:lnbc(outfile))//'.bin',
     & access='direct',status='unknown',form='unformatted',recl=4*nsig)
      irec=0
      do ilev=1,nlev
         read(lun_mav,*) (vv(icol),icol=1,ncol)
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
c   AbsW is the difference: Full LM  -  1'st order abscoeffs (CO2)
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
         endif

         write(lun_wbin,rec=ilev) (sngl(AbsY(i)),i=1,nsig)
         do i=1,nsig
            sig=sigmin+(i-1)*dsig
         write(lun_vac,123)ilev,vv(2),vv(3),sig,AbsV(i),AbsY(i),AbsW(i)
 123     format(i4,f8.2,f8.5,f12.5,3(1pe12.3))
         end do  ! i=1,nsig
         irec=irec+nsig
      end do   ! ilev=1,nlev
      close(lun_vac)
      go to 1   !  Loop to next "Next Spectrum"
99    close(lun_mav)
      close(lun_wbin)
      stop
      end
