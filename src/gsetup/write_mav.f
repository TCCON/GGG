      subroutine write_mav(z,t,p,d,vmr,nlev,lunw_mav,lncell,
     & t_cell,p_cell,vmr_cell,gas_in_cell,vmrlabel,isofile)
c
c  Writes a single set of Z/T/P/VMR profiles to LUNW_MAV.
c  A complete .mav file typically contains multiple sets of profiles.
c
c  Does NOT write the header info to LUNW_MAV.
c  This must be done from the calling program.
c
c  Inputs:
c    Z(NLEV)        R*4   Altitudes (km) of the levels
c    T(NLEV)        R*4   Temperatures (K) of the levels
c    P(NLEV)        R*4   Pressures (atm) of the levels
c    D(NLEV)        R*4   Densities (cm-3) of the levels
c  VMR(NGAS,NLEV)   R*4   DMFs of atmospheric gases
c    NLEV           I*4   Number of Levels
c    LUNW_MAV       R*4   Logical Unit Number of .mav file
c    LNCELL         I*4   Number of Cells
c    T_Cell(NCELL)  R*4   Temperatures of cells
c    P_Cell(NCELL)  R*4   Pressures of cells
c  VMR_Cell(NCELL)  R*4   WMFs of gases in cells
c  GAS_Cell(NCELL)  R*4   Gas ID of cells
c     VMRLABEL      C**   Column label string from .vmr file
c     ISOFILE       C**   Path to isotopologs.dat file
c
c Outputs:
c   .mav file containing T, P, WMFs
c
c  Isotopic fractionation is imparted to the VMR profiles.
c  The observed isotopic ratio R(z) at altitude z is related to
c  the isotopic ratio in the troposphere, R(0), by the equation
c     R(z) = R(0).f(z)^(lnfrd/1000)
c  see for example Griffith et al., GRL, 27, 2485-2488, [2000]
c  f(z) is the fraction of gas remaining at z: f=VMRp(z)/VMRp(0)
c  lnfrd is the relative enrichment factor (per mil).
c
c  Since f < 1 in the stratosphere, a negative value of lmfrd
c  causes the frationation to increase with altitude. For example,
c     0.9^(50/1000) = 0.995;  0.90^(-50/1000) = 1.005
c
c  Let VMRp(z) is the mole fraction of the parent isotopolog,
c  Let VMRi(z) is the mole fraction of isotopolog i (divided by IUPAC ratio).
c  Since R(z)=VMRi(z)/VMRp(z)
c     VMRi(z)/VMRp(z) = [VMRi(0)/VMRp(0)]*[VMRp(z)/VMRp(0)]^(lnfrd/1000)
c  which simplifies to
c     VMRi(z) = VMRi(0)*[VMRp(z)/VMRp(0)]^(1+lnfrd/1000)
c  Since VMRi(0)=(1+delta/1000)*VMRp(0)
c     VMRi(z) = (1+delta/1000)*VMRp(0)*[VMRp(z)/VMRp(0)]^(1+lnfrd/1000)
c     VMRi(z) = (1+delta/1000)*VMRp(z)^(1+lnfrd/1000) / VMRp(0)^(lnfrd/1000)
c
c  Isotopic fractionation of H2O is treated as a special case.
c  In gsetup.f, HDO is set to H2O*0.15*(8.0+log10(h2ovmr(ilev)))
c  In this subroutine, the 2h2o and 3h2o isotoplog fractionations
c  are set to 1/8 and 1/16 of that of HDO, respectively.
c
c  GASINDEX maps the species (in isotopolog.dat) into gases
c    kgas=gasindex(jspeci)
c  e.g.  gasindex(56)=14  ! HF
c  e.g.  gasindex(57)=14  ! DF
c  e.g.  gasindex(58)=15  ! HCl
c
c  SAGINDEX maps the gases into species (first isotopolog)
c    jspeci=sagindex(jvmr)
c  e.g.  sagindex(14)=56  ! HF
c  e.g.  sagindex(15)=58  ! HCl
c
c  Hence GASINDEX and SAGINDEX are the inverse of each other.
c    gasindex(sagindex(jvmr))=jvmr

      implicit none

      real*4  zero,zdiff
      parameter(zero=0.0)

      include "../gfit/ggg_int_params.f"

      integer*4 lunw_mav,nspeci,jspeci,nlev,i,j,kspeci,
     & i10,
     & istat,jvmr,kvmr,nvmr,lunr_iso,
     & mcell,jcell,lncell,gas_in_cell(lncell)
      parameter (kspeci=230,kvmr=80,lunr_iso=57,mcell=8)

      integer*4 gasindex(kspeci),sagindex(kvmr),
     & nlhead,ncol_iso,ip_h2o,
     & kgas,kiso,dumint,ifail,idum,icode,
     & dumintarr(mvmode),nmode,
     & speci_cell(mcell),nspeci_cell(mcell)

      real*4 z(nlev+1),t(nlev+1),p(nlev+1),d(nlev+1),vmr(mgas,nlev+1),
     &  delta(kspeci),lnfrd(kspeci),fia,dumreal,
     &  dumrealarr(mvmode),boltzmann,
     &  vmrh2o1,vmrh2o2,vmrh2o3,
     &  xmrh2o1,xmrh2o2,xmrh2o3,
     &  t_cell(lncell),p_cell(lncell),vmr_cell(lncell)

      character isofile*(*),
     &     shortname(kspeci)*(8), fullname(kspeci)*(11),
     &     vmrlabel*(*),vmrgases(kvmr)*10,speci_id*20,iso_fmt*80

      ip_h2o=1  ! column mumber of H2O in VMR array
      boltzmann=1.38066E-23
      data sagindex/kvmr*0/

      idum=mauxcol   ! Prevent compiler warning (unused variable)
      idum=mcolvav   ! Prevent compiler warning (unused variable)
      idum=mfilepath ! Prevent compiler warning (unused variable)
      idum=mlev      ! Prevent compiler warning (unused variable)
      idum=mrow_qc   ! Prevent compiler warning (unused variable)
      idum=mspeci    ! Prevent compiler warning (unused variable)
      idum=ncell     ! Prevent compiler warning (unused variable)
      idum=nchar     ! Prevent compiler warning (unused variable)

      if(lncell.gt.mcell) stop 'write_mav: lncell > mcell'
      open(19,file='write_mav.rpt')
      write(19,*)'write_mav.rpt (force over-write of previous version)'
      do jcell=1,lncell
         speci_cell(jcell)=0
         nspeci_cell(jcell)=0
c     Index gas names in vmrlabel to gas names in "isotopologs.dat"
c      write(6,*)' Opening ',isofile
         open(lunr_iso,file=isofile,status='old')
         read(lunr_iso,*) nlhead,ncol_iso
         read(lunr_iso,'(7x,a)') iso_fmt
         do j=3,nlhead
            read(lunr_iso,*) 
         end do
         do jspeci=1,kspeci
            call read_isotopolog(lunr_iso,iso_fmt,kgas,kiso,
     &      shortname(jspeci),speci_id,icode,fia,delta(jspeci),
     &      lnfrd(jspeci),dumint,dumreal,dumreal,dumreal,dumrealarr,
     &      dumintarr,nmode,mvmode,istat)
            if(istat.ne.0) goto 77
            if(kgas.eq.gas_in_cell(jcell)) then
               speci_cell(jcell)=jspeci
               nspeci_cell(jcell)=nspeci_cell(jcell)+1
            endif
            call lowercase(shortname(jspeci))
            write(fullname(jspeci),'(i2,a)') kiso,shortname(jspeci)
            gasindex(jspeci)=0
c            write(*,*) jspeci,kgas,shortname(jspeci), speci_id,
c     &           fia,delta(jspeci),lnfrd(jspeci),
c     &           dumint, dumreal, nmode
         end do
         write(6,*) 'write_mav: Increase parameter KSPECI'
         stop
 77      close(lunr_iso)
         nspeci=jspeci-1
      end do       !  jcell=1,lncell

      call lowercase(vmrlabel)
      call substr(vmrlabel,vmrgases,kvmr,nvmr)
      nvmr=nvmr-1               !  The first column of the vmr file is Altitude.
c      write(*,*)'nspeci,nvmr=',nspeci,nvmr
      if(nvmr.gt.mgas) stop ' GSETUP: NVMR > MGAS'
c     
c  Check whether all the gas names in "isotopologs.dat" were also
c  found in the .vmr file.
      ifail=0
      call clistindex(nspeci,fullname,nvmr,vmrgases(2),gasindex)
      call clistindex(nspeci,shortname,nvmr,vmrgases(2),gasindex)
      do jspeci=1,nspeci
         if(gasindex(jspeci).le.0) then
            write(19,*) ' write_mav: Warning: '//
     &      shortname(jspeci)//' not found in .vmr file '
            ifail=1
         endif
c         write(*,*)ifail,jspeci,'gasindex=',gasindex(jspeci)
      end do
c
c  Check whether all the gas names in the .vmr file were also
c  found in "isotopologs.dat"
      call clistindex(nvmr,vmrgases(2),nspeci,fullname,sagindex)
      call clistindex(nvmr,vmrgases(2),nspeci,shortname,sagindex)
      do jvmr=1,nvmr
c       write(*,*)'sagindex=',jvmr,sagindex(jvmr),gasindex(sagindex(jvmr))
         if(sagindex(jvmr).le.0) write(19,*) ' write_mav: Warning: '//
     &   vmrgases(1+jvmr)//' not found in isotopologs.dat file: '
      end do
      close(19)
      if(ifail.eq.1) then
         write(*,*)' Error: Gases missing from .vmr file'
         write(*,*)' See write_mav.rpt for details'
         stop 'GSETUP Error'
      endif
c     
c   Convert atmospheric DMFs into WMFs in place
      call dmf2wmf(mgas,mgas,nlev,ip_h2o,vmr)

c     Output model information (SUNRUN.MAV)
c      write(lunw_mav,'(3i4)') 2, nspeci+4, nlev+lncell
      write(lunw_mav,'(a,230a12)')
     & ' Height  Temp     Pres       Density    ',
     &   (fullname(jspeci),jspeci=1,nspeci)

c  Find level (110) closest to 10 km (tropopause).
      zdiff=z(1)-10.
      i10=1
      do i=2,nlev
         if(abs(z(i)-10.0).lt.abs(zdiff)) then
            zdiff=z(i)-10.0
            i10=i
         endif
      enddo
c      write(*,*) 'i10 = ',i10,' z(i10) = ',z(i10)
c
c First NCELL lines of .mav file contains cell information
      do jcell=1,lncell
c         write(*,*)jcell,speci_cell(jcell),nspeci_cell(jcell)
         write(lunw_mav,'(2f7.2,2(1pe12.4),230e12.4)')
     $   -11+1.1*float(jcell),t_cell(jcell),p_cell(jcell)/1013.25,
     &   0.0001*p_cell(jcell)/boltzmann/t_cell(jcell),
     &   (zero,j=1,speci_cell(jcell)-nspeci_cell(jcell)),
     &   (vmr_cell(jcell),j=1,nspeci_cell(jcell)),
     &   (zero,j=speci_cell(jcell)+1,nspeci)
      end do      ! jcell=1,lncell

c Last NLEV lines of .mav file contain atmospheric info
c     Write tropospheric fractions
      do i=1,i10
         vmrh2o1=vmr(1,i)
         vmrh2o2=(vmr(1,i)**0.875)*(vmr(49,i)**0.125)
         vmrh2o3=(vmr(1,i)**0.9375)*(vmr(49,i)**0.0625)
c         vmr(6,i)=45E-09  ! Mars Korablev
c         write(86,*) z(i),t(i),p(i),d(i),vmrh2o1,vmr(2,i)
         write(lunw_mav,'(2f7.2,2(1pe12.4),230e12.4)')
     &    z(i),t(i),p(i),d(i),vmrh2o1,vmrh2o2,vmrh2o3,
     &   (vmr(gasindex(j),i)*(1.0+0.001*delta(j)),j=4,nspeci)
      end do      ! i=1,i10
c  Write stratospheric enrichments.
c  Below, vmrh2o is the vmr at i10, whereas xmvh2o is the vmr
c  at the current altitude.
      do i=i10+1,nlev
         xmrh2o1=vmr(1,i)
         xmrh2o2=(vmr(1,i)**0.875)*(vmr(49,i)**0.125)
         xmrh2o3=(vmr(1,i)**0.9375)*(vmr(49,i)**0.0625)
c         vmr(6,i)=45E-09  ! Mars Korablev
         write(lunw_mav,'(2f7.2,2(1pe12.4),230e12.4)')
     &   z(i),t(i),p(i),d(i),
     &   (1.0+0.001*delta(1))*xmrh2o1**(1.0+0.001*lnfrd(1))/
     &   vmrh2o1**(0.001*lnfrd(1)),
     &   (1.0+0.001*delta(2))*xmrh2o2**(1.0+0.001*lnfrd(2))/
     &   vmrh2o2**(0.001*lnfrd(2)),
     &   (1.0+0.001*delta(3))*xmrh2o3**(1.0+0.001*lnfrd(3))/
     &   vmrh2o3**(0.001*lnfrd(3)),
     &   ((1.0+0.001*delta(j))*vmr(gasindex(j),i)**
     &   (1.0+0.001*lnfrd(j))/vmr(gasindex(j),i10)**
     &   (0.001*lnfrd(j)),j=4,nspeci)
      end do      ! i=i10+1,nlev
      return
      end
