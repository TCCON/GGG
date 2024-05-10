c  write_aux.f
c  Writes out a nice version of the mav file for data users. 

      implicit none
      include "../gfit/ggg_int_params.f"

      integer*4 lun_mav,lun_out,nlev,nspeci,lnbc,lcolon,ldot,i,idum,
     & j,ii,lspace,ii1,jj,ico2,ico,ich4,in2o,ih2o,ihdo,ihf,io2,flag,
     & nlhead
      parameter (lun_mav=12)      ! input file (.mav)
      parameter (lun_out=14)      ! output file (.map)
c      parameter (mvmr=28000)
      real*4 z(mlev),t(mlev),p(mlev),d(mlev),vmr(mspeci*mlev),
     & gravity,oblat,wmf_h2o
      character version*64,mavfile*80,outfile*80,mavstring*64,
     & runlabmav*57,string*48,head*1800,hdr(mlev)*(11),string1*28,
     & mav_version*64,vmrfile*1800,modfile*1800

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
     &' WRITE_AUX                Version 1.13      2019-02-13      GCT '
      write(*,*) version

      oblat=0.0d0
      if (iargc() == 0 ) then
         write(*,'(a)') 'Enter name of .mav file'
         read(*,'(a)')mavfile
      elseif (iargc() == 1) then
         call getarg(1, mavfile)
      else
c         stop 'Usage: $gggpath/bin/write_aux mavfilename latitude'
         stop 'Usage: $gggpath/bin/write_aux mavfilename'
      endif

c      if (iargc() == 0) then
c         write(*,'(a)')
c     &    'Enter latitude of site (default: 0 degrees)'
c         read(*,'(f8.4)')oblat
c      endif

c  Read the entire contents of the .mav file
      open(lun_mav,file=mavfile,status='old')
      read(lun_mav,'(a)')mav_version
      read(lun_mav,'(14x,a)')string
      lcolon=index(string,':')
      read(string(lcolon+1:),'(a)')runlabmav

      do i=1,9999
         call read_mav_head(lun_mav,nlev,nspeci,nlhead)
         call read_mav_aux(lun_mav,nlhead,nlev,nspeci,z,t,p,d,vmr,
     & oblat,vmrfile,modfile,head)
c        write(*,*)modfile(:lnbc(modfile))
         ldot=index(modfile,'.')
         outfile=runlabmav(1:2)//modfile(ldot-20:ldot-10)//'.map'
c The outfile will be in the form: pa20040721.map
c         write(*,*)outfile
c         write(*,*)runlabmav
c         write(*,*)head
c         write(*,*)mav_version
c         outfile=runlabmav(:lr-13)//'.map'
c         lr=lnbc(runlabmav)

c  Figure out where the date ends (usually has an 's' for TCCON spectrum names).
c         lr=index(runlabmav,'s')
c         if (lr.lt.3) lr=11 ! To avoid problems with station names with s in the site id or no site id at all
c         write(*,*)lr
c         outfile=runlabmav(:lr-1)//'.map'
c         write(*,*)'old=',outfile

c  Figure out what is in the mav file and select out the appropriate columns
         ii1=1
         jj=1
         ico2=0
         ico =0
         ich4=0
         in2o=0
         ih2o=0
         ihdo=0
         ihf =0
         io2=0
         flag=0
         do ii=1,99999
c            write(*,*)ii1
            if(head(ii1:ii1).eq.' ') then
               ii1=ii1+1
            else
               lspace=index(head(ii1:),' ')+ii1
               if(index(head(ii1:),' ')+ii1.ge.lnbc(head)) then
                  hdr(jj)=head(ii1:)
                  flag=1
               else
                  hdr(jj)=head(ii1:lspace)
               endif
               if(hdr(jj).eq.'1co2') then
                  ico2 = jj
               elseif(hdr(jj).eq.'1co') then
                  ico  = jj
               elseif(hdr(jj).eq.'1h2o') then
                  ih2o = jj
               elseif(hdr(jj).eq.'1hdo') then
                  ihdo = jj
               elseif(hdr(jj).eq.'1n2o') then
                  in2o = jj
               elseif(hdr(jj).eq.'1ch4') then
                  ich4 = jj
               elseif(hdr(jj).eq.'1hf') then
                  ihf  = jj
               elseif(hdr(jj).eq.'1o2') then
                  io2  = jj
               endif
               ii1=lspace+1
               jj=jj+1
            endif
            if(flag.eq.1) then
c               write(*,*)'flag is 1'
               goto 77
            endif
         enddo ! do ii=1,99999 parse headers
77       continue

c         write(*,*)outfile
c         write(*,*)'ico2,ico,ih2o,ihdo,ich4,in2o=',
c    &               ico2,ico,ih2o,ihdo,ich4,in2o
         if(ico2.eq.0) then
            write(*,*)'Warning, 1co2 was not found.'
         endif
         if(ico .eq.0) then
            write(*,*)'Warning, 1co was not found.'
         endif
         if(ich4.eq.0) then
            write(*,*)'Warning, 1ch4 was not found.'
         endif
         if(in2o.eq.0) then
            write(*,*)'Warning, 1n2o was not found.'
         endif
         if(ih2o.eq.0) then
            write(*,*)'Warning, 1h2o was not found.'
         endif
         if(ihdo.eq.0) then
            write(*,*)'Warning, 1hdo was not found.'
         endif
         if(ihf .eq.0) then
            write(*,*)'Warning, 1hf was not found.'
         endif
         if(io2 .eq.0) then
            write(*,*)'Warning, 1o2 was not found.'
         endif

c Determine which columns to output
         string1=hdr(ih2o)(2:4)//','//hdr(ihdo)(2:4)//','
     &   //hdr(ico2)(2:4)//','//hdr(in2o)(2:4)//','//hdr(ico)(2:3)//','
     &   //hdr(ich4)(2:4)//','//hdr(ihf)(2:3)//','//hdr(io2)(2:3)
     
c         write(*,*)string1

c Write new .map file
         open(lun_out,file=outfile,status='unknown')
         write(lun_out,'(1x,i2,1x,i2)')11,12
         write(lun_out,'(1x,a)')outfile(:lnbc(outfile))
         write(lun_out,'(a)')mav_version
         write(lun_out,'(a)')version
         write(lun_out,'(1x,a)')
     & 'Please see https://tccon-wiki.caltech.edu for a complete'//
     & ' description of this file and its usage.'
         write(lun_out,'(1x,a)')
     & 'Avogadro (molecules/mole): 6.0221415E+23'
         write(lun_out,'(1x,a)')
     & 'Mass_Dry_Air (kg/mole): 28.9644E-03'
         write(lun_out,'(1x,a)')
     & 'Mass_H2O (kg/mole): 18.01534E-03'
         write(lun_out,'(1x,a,f7.3)')
     & 'Latitude (degrees): ',oblat
         write(lun_out,'(1x,a,a,a)')'Height,Temp,Pressure,Density,',
     & string1(:lnbc(string1)),',gravity'

c Match units to string1
         write(lun_out,'(1x,a,a)')'km,K,hPa,molecules_cm3,',
     & 'parts,parts,ppm,ppb,ppb,ppb,ppt,parts,m_s2'

c    Note that the pe11.3 format without 0p affects the subsequent f11.3 format.
c    It multiplies it by 10!
10       format ((2(f7.2,","),1p,4(e10.3,","),0p,1(f7.2,","),1p,
     & (e10.3,","),0p,(f9.3,","),(f7.1,","),(f8.2,","),(f7.4,","),f6.3))
         do j=1,nlev
            if(z(j).ge.0) then ! Skip cell levels
               wmf_h2o = 1 - vmr(ih2o-4+(j-1)*nspeci)/0.997317 ! Approximate the wet mole fraction of H2O
               write(lun_out,10)z(j),t(j),
     &         p(j)*1013.25,d(j),
     &         vmr(ih2o-4+(j-1)*nspeci),                     ! H2O WMF in parts
     &         vmr(ihdo-4+(j-1)*nspeci),                     ! HDO WMF in parts
     &         vmr(ico2-4+(j-1)*nspeci)*1e6,                 ! CO2 WMF in ppm
     &         vmr(in2o-4+(j-1)*nspeci)*1e9,                 ! N2O WMF in ppb
     &         vmr(ico -4+(j-1)*nspeci)*1e9,                 ! CO  WMF in ppb
     &         vmr(ich4-4+(j-1)*nspeci)*1e9,                 ! CH4 WMF in ppb
     &         vmr(ihf -4+(j-1)*nspeci)*1e12,                ! HF  WMF in ppt
     &         vmr(io2 -4+(j-1)*nspeci),                     ! O2  WMF in parts
c     &         (vmr(ih2o-4+(j-1)*nspeci)/wmf_h2o),           ! H2O DMF in parts
c     &         (vmr(ihdo-4+(j-1)*nspeci)/wmf_h2o),           ! HDO DMF in parts
c     &         (vmr(ico2-4+(j-1)*nspeci)/wmf_h2o)*1e6,       ! CO2 DMF in ppm
c     &         (vmr(in2o-4+(j-1)*nspeci)/wmf_h2o)*1e9,       ! N2O DMF in ppb
c     &         (vmr(ico -4+(j-1)*nspeci)/wmf_h2o)*1e9,       ! CO  DMF in ppb
c     &         (vmr(ich4-4+(j-1)*nspeci)/wmf_h2o)*1e9,       ! CH4 DMF in ppb
c     &         (vmr(ihf -4+(j-1)*nspeci)/wmf_h2o)*1e12,      ! HF  DMF in ppt
c     &         (vmr(io2 -4+(j-1)*nspeci)/wmf_h2o),           ! O2  DMF in parts
     &         gravity(oblat,z(j))                           ! gravity in m/s^2
            endif
         enddo ! do j=1,nlev
         close(lun_out)

         read(lun_mav,'(a)',end=66) mavstring
         if(mavstring(1:14).eq.'Next Spectrum:') then
            read(mavstring(15:),'(a)') runlabmav
         else
            write(*,*) mavstring
            write(*,*)'Failed to find Next Spectrum string'
            stop
         endif ! Next Spectrum
      enddo ! do 1,9999
66    continue
      close(lun_mav)
      end
