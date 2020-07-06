c   Provides user interface to create stand-alone apriori vmr profiles
c
c  Inputs:
c     reference.vmr      ! Initial set of vmr profiles
c     saap_input.lst     ! List of dates/latitudes for which profiles are needed
c
c  Outputs:
c     saap_YYYYMMDD_LL_ZZZZZ.out   !  VMR profiles
c
c   YYYY = Year
c   MM   = Month
c   DD   = Date
c   LL   = Latitude (deg)
c   ZZZZZ= Tropopause altitude (m)

      implicit none

      integer*4
     & nlheader, ncol,
     & nlev,ilev,
     & mgas,
     & idoy,iyr,
     & lunr_lst,
     & lunr_vmr,
     & lunw_sum,
     & jgas,
     & ngas,
     & imm,idd,jul,jul0 
      parameter (lunr_lst=69)    ! for reading the .lst file
      parameter (lunr_vmr=73)    ! for reading the .vmr file (readvmrFC)
      parameter (lunw_sum=83)    ! for writing the .out file (readvmrFC)
      parameter (nlev=71)        ! maximum number of atmospheric levels
      parameter (mgas=80)        ! maximum number of gases

      real*4 
     & gradlat(mgas),       ! Latitude Gradients
     e strend(mgas),        ! Secular Trends
     & seacycle(mgas),      ! Seasonal Cycle
     & refvmr(mgas,nlev),   ! buffer for reference vmr's
     & apvmr(mgas,nlev),    ! buffer for a priori vmr's
     & z(nlev)              ! altitudes of levels (km)

      real*8 fryr,reflat_vmr,date_mod,date_vmr,ztrop_vmr

      real*8 ztrop_gct,
     & oblat                ! observation latitude (deg).
        
      character
     & cc*1,
     & vmrlabel*1024,       ! column labels from vmr file
     & outfile*120,
     & saap_lst*80,         ! file name for list of locations
     & version*64           !current program version

      version=
     & ' stand_alone_a_priori     Version 0.09     2016-10-23    GCT '

c  choose an observation geometry
      write(6,*) version
c------------------------------------------------------------------
      do ilev=1,nlev
         z(ilev)=float(ilev-1)
      end do
c---------------------------------------------
c  Open the file which will contain the slant paths
c  Get the header/runlog information pertaining to RUNLAB
      if (iargc() == 0 ) then
         write(*,'(a)') 'Enter name of .lst file (e.g. saap_input.lst)'
         read(*,'(a)')saap_lst
      elseif (iargc() == 1) then
         call getarg(1, saap_lst)
      else
         stop 'Usage: $gggpath/bin/stand_alone_a_priori saap_input.lst'
      endif

c      open(lunr_lst,file='saap_input.lst',status='old')
      open(lunr_lst,file=saap_lst,status='old')
      read(lunr_lst,*) nlheader,ncol
      read(lunr_lst,*)
      do  !  Loop over reading .lst file
         read(lunr_lst,*,end=99) iyr,imm,idd,oblat,ztrop_gct
         write(*,*)  iyr,imm,idd,oblat,ztrop_gct
         ztrop_gct=ztrop_gct/1000   ! convert m to km
         if(oblat.gt.0.0) then
            cc='N'
         else
            cc='S'
         endif
         call julian(iyr,imm,idd,jul)
         call julian(iyr,1,1,jul0)
         idoy=jul-jul0+1
         date_mod=iyr+idoy/365.25 !from the runlog

         call read_refvmrs(lunr_vmr,'summer_35N.vmr',nlev,z,mgas,'mmmm',
     &   vmrlabel,refvmr,ngas,strend,gradlat,seacycle,
     &   reflat_vmr,date_vmr,ztrop_vmr)

c If there's only 1 level (i.e. lab), or the vmr file for that particular
c measurement already exists, don't try to modify the .vmr file
         call resample_vmrs_at_effective_altitudes(nlev,z,mgas,
     &   ngas,refvmr,ztrop_gct,ztrop_vmr,oblat,reflat_vmr,apvmr)
         call apply_vmr_latitude_gradients(nlev,z,mgas,ngas,gradlat,
     &   apvmr,ztrop_gct,reflat_vmr,oblat,apvmr)
         call apply_secular_trends(nlev,z,mgas,ngas,strend,apvmr,
     &   ztrop_gct,reflat_vmr,oblat,date_mod,date_vmr,apvmr)
         fryr=date_mod-int(date_mod)
         call apply_seasonal_cycle(nlev,z,mgas,ngas,seacycle,
     &   ztrop_gct,oblat,fryr,apvmr)
c         do ilev=1,nlev
c            z8=dble(z(ilev))
c            do jgas=1,ngas
c               apvmr(jgas,ilev)=apvmr(jgas,ilev)*
c     &         compute_seasonal_cycle(jgas,z8,ztrop_gct,oblat,fryr)
c            end do
c         end do

c  Output vmr info (no isotopic fractionation)
         write(outfile,'(a5,i4,2i2.2,a1,i2.2,a1,a1,i5.5,a)')'saap_',
     &   iyr,imm,idd,'_',nint(abs(oblat)),cc,'_',nint(1000*ztrop_gct),
     &   '.out'
         open(lunw_sum,file=outfile,status='unknown')
         write(lunw_sum,*)2,1+ngas
         write(lunw_sum,'(a)') 
     & '   Z    H2O      CO2       O3        N2O       CO        CH4  
     &     O2    NO         SO2        NO2        NH3        HNO3  
     &     OH         HF      HCl        HBr        HI         ClO  
     &     OCS        H2CO       HOCl     HO2        H2O2       HONO 
     &     HO2NO2     N2O5       ClNO3      HCN     CH3F       CH3Cl 
     &     CF4        CCl2F2     CCl3F      CH3CCl3    CCl4    CCOF2 
     &     COFCl      C2H6       C2H4       C2H2       N2    CCHClF2
     &     COCl2      CH3Br      CH3I       HCOOH      H2S    CCHCl2F
     &     HDO        SF6        F113       ClCN       F142b  Cdust_m
     &     PH3        CH3OH      CH3SH      CH3CHO     CH3CN      PAN
     &     NF3        ClOOCl     ClClO2     ClOClO     CHF3     f141b
     &     CCH3COOH    cirrus6    cirrus15   C3H8    D2O    SAV   Air'
         do ilev=1,nlev
            write(lunw_sum,'(f6.1,80e12.5)') z(ilev),
     &      (apvmr(jgas,ilev),jgas=1,ngas)
         end do
         close(lunw_sum)
c------------------------------------------------------------------
      end do     !  Loop over reading .lst file
 99   close(lunr_lst)
      stop
      end

