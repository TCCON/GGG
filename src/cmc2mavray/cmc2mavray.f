c  Reads the .cmc files 
c  From the information therein
c  produces .mav and .ray files.
c  CMC = Cell Measurement Conditions
c  Used in the analysis of lab spectra.
c
c  cmc2mavray   V 1.06   2016-09-04   GCT'
c  Since HITRAN multiplies all line intensities by the
c  fractional isotopic abundances, even for the parent
c  isotopolog, then all the vmr values written to the
c  .mav file for use by GFIT represent the total (all
c  isotoplogs) gas amount, even though they are 
c  attributed to individual isotopologs. For example,
c  the vmr attributed to 12CH4 and 13CH4 are identical,
c  even though the latter is 90x smaller in reality.
c
c  This is convenient when gas samples are unfractionated
c  since it allows calculations to be performed in blissful
c  ignorance of the isotopic abundances. The fact that
c  the 13CH4 VMRs are 90x too large is compensated by the
c  line strengths being 90x too small.
c  But when the gas sample is isotopically enriched, the
c  isotopic ratios built into HITRAN must first be removed
c  before the true vmrs of the isotopologs can be used.
c  
c  There are two different ways that cmc2mavray could
c  handle this.
c  1) Input the true partial pressure of each isotopolog
c   into the .cmc file.  cmc2mavray would then divide
c   them by the isotopic abundances built into HITRAN
c   in order to undo the intensity fudging.
c  2) Input the partial pressures already divided by
c   the isotopic fractions into the .cmc file. This
c   means that in highly enriched cases, the partial
c   pressures can exceed the total pressure.
c
c  We implemented Option 1 despite the disadvantage of
c  needing an additional input file (isotopologs.dat)
c  containing the isotopic fractionations that were
c  applied to the line intensities.  Option 1 has the
c  advantage of allowing meaningful checks that the
c  isotopic partial pressures don't exceed the total.

      implicit none
      integer*4 j,lc,lh,jcell,ncell,mcell,lrt,
     & kgas_iso,kiso_iso,icode,istat,lnblnk,
     & mspeci_iso,nspeci_iso,jspeci,igas,mgas,
     & ifirst,ilast,bytepw,possp,iyr,iset,
     & nlhead,ncol,
     & lunr_iso,lunr_cmc,lunw_ray,lunw_mav,lunr_rlg,irec
      parameter (mgas=16,mcell=5,mspeci_iso=200,
     & lunr_iso=11,lunr_cmc=12,lunw_ray=14,lunw_mav=15,lunr_rlg=16)
      integer*4 ngas(mcell),kgas_cmc(mgas,mcell),
     & kiso_cmc(mgas,mcell),
     & gasindex(mspeci_iso),
     & specindex(80)
      real*4 pp(mgas,mcell),cleng(mcell),ctemp(mcell),tpp,
     & cpres(mcell),dens,vmr_mav(mspeci_iso),zero,fia(mspeci_iso)
      character specname_cmc*57,cmcfile*256,specname_rlg*57,
     & data_fmt_read_rl*256, col_labels_rl*320,apf*3,
     & gggdir*256,dl*1,iso_fmt*80,speci_id*20,col1*1,
     & format_ray*80,version*48,gasname*8,header*2050

      real*8 graw,
     & oblat,            ! observation latitude (deg).
     & oblon,            ! observation longitude (deg).
     & obalt,            ! observation altitude (km)
     & zpdtim,           ! Time of ZPD (UT hours)
     & asza,             ! astronomical solar zenith angle (unrefracted)
     & zenoff,           ! zenith pointing offset
     & azim,             ! azimuth angle
     & osds,             ! Observer-Sun Doppler Stretch (ppm)
     & fovi,             ! Internal angular diameter of FOV (radians)
     & fovo,             ! External angular diameter of FOV (radians)
     & amal,             ! angular misalignment of interferometer
     & zoff,             ! zero-level offset
     & snr,              ! Signal-to-Noise Ratio
     & tins,             ! Temperature INSide the instrument
     & pins,             ! Pressure INSide the instrument
     & hins,             ! Humidity INSide the instrument
     & tout,             ! Temperature OUTside the instrument
     & pout,             ! Pressure OUTside the instrument
     & hout,             ! Humidity OUTside the instrument
     & sia,              ! Solar Intensity Average (arbitrary units)
     & fvsi,             ! Fractional Variation in Solar Intensity
     & wspd,             ! Wind Speed
     & wdir,             ! Wind Direction
     & aipl,             ! Airmass-Independent Path Length (km)
     & lasf,             ! laser frequency (e.g. 15798.03 cm-1)
     & wavtkr,           ! suntracker operating frequency (9900 cm-1)
     & opd               ! Optical path difference (cm) of interferogram


      version=' cmc2mavray    Version 1.19    2020-01-26   GCT '

c  Environment specification: 
      call get_ggg_environment(gggdir, dl)
      lrt=lnblnk(gggdir)     !Length of gggdir
c
      zero=0.0
      write(format_ray,'(a,i2.2,a)') '(a57,',7+mcell,'(1x,f13.7))'

c  Get name of input file 
      lc=0
      do while (lc .le. 0)
         if (iargc() == 0) then
            write(6,'(a)') 'Input File (e.g. h2o_kp.cmc): '
            read(*,'(a)')cmcfile
         elseif (iargc() == 1) then
            call getarg(1, cmcfile)
         else
            stop 'Use: $gggpath/bin/cmc2mavray cmcfile '
         endif
         lc=lnblnk(cmcfile)
      end do   !  while (lc .le. 0)

      header='  Height  Temp   Pres   Density      '
      lh=lnblnk(header)+5
c  Read in names of isotopomers
      open(lunr_iso,file=gggdir(:lrt)//'isotopologs/isotopologs.dat',
     & status='old')
      read(lunr_iso,*) nlhead,ncol
      read(lunr_iso,'(a)') iso_fmt
      do j=3,nlhead
         read(lunr_iso,*)               ! Skip header info + column labels
      end do
      do jspeci=1,mspeci_iso
1        read(lunr_iso,iso_fmt,end=77)col1,kgas_iso,kiso_iso,
     &   gasname,speci_id,icode,fia(jspeci)
c         write(*,*) jspeci,kgas_iso,kiso_iso,gasname
         if(col1.eq.':' .or. col1.eq.';') go to 1

c         write(header(lh+1:lh+11),'(i3,a8)') mod(kiso_iso,10),gasname
         write(header(lh+1:lh+11),'(i3,a8)') kiso_iso,gasname
         lh=lh+11
         specindex(kgas_iso)=jspeci-kiso_iso
         gasindex(jspeci)=kgas_iso
         if( specindex(kgas_iso).lt.0) then
            write(*,*) jspeci,kgas_iso,kiso_iso, specindex(kgas_iso)
            stop 'specindex < 0'
         endif
      end do
      read(lunr_iso,*,end=77)
      stop ' GFIT: Number of species in ISOTOPOLOG.DAT exceeds MSPECI'
77    close(lunr_iso)
      nspeci_iso=jspeci-1
      write(*,*)'nspeci_iso=',nspeci_iso
c      write(*,*)'gasindex(179)=',gasindex(179)
c      write(*,*)'specindex(76)=',specindex(76)
      call lowercase(header)
c      write(*,*) header
c      write(*,*) specindex

      open(lunw_ray,file=cmcfile(:lc-4)//'.ray',status='unknown')
      write(lunw_ray,*) 4,7+mcell
      write(lunw_ray,'(a)') version
      write(lunw_ray,'(a)') 'format='//format_ray
      write(lunw_ray,'(a,5(a6,i1.1))')
     & ' SpectrumName  Zobs        Pobs        ASZA         Bend        
     & FOV        Zmin     ',('  Cell',j,j=1,mcell)

      open(lunw_mav,file=cmcfile(:lc-4)//'.mav',status='unknown')
      write(lunw_mav,'(a)') version

c  Open the runlog to check for a one-to-one correspondance
c  between spectra in runlog and cmc file.
      open(lunr_rlg,file=gggdir(:lrt)//'runlogs/lab/'//
     & cmcfile(:lc-4)//'.lrl',status='old')
      call read_runlog_header(lunr_rlg,data_fmt_read_rl,col_labels_rl)

      open(lunr_cmc,file=cmcfile,status='old')

      write(66,*)4,7
      write(66,*)'  Spectrum_Name          Length    Temp  Tot_Pres   PP
     &      VMR    v_start  v_end'
      write(66,*)'                          (m)      (K)    (Torr)   (To
     &rr)  parts (cm-1)  (cm-1)'
      write(66,*)'  Spectrum_Name          Length    Temp  Tot_Pres   PP
     &      VMR    v_start  v_end'
      write(*,*)'  Spectrum_Name          Length    Temp  Tot_Pres    PP
     &      VMR    v_start  v_end'
      do irec=1,999999
22       call read_cmc(lunr_cmc,specname_cmc,mcell,ncell,
     &   kgas_cmc,kiso_cmc,cleng,ctemp,cpres,pp,mgas,ngas,istat)
         if(istat.eq.3) go to 22  ! colon in col 1
c         write(*,*) 'cpres = ',(cpres(j),j=1,ncell)
         if(istat.ne.0) go to 99
c
         call read_runlog_data_record(lunr_rlg,data_fmt_read_rl,
     &   col1,specname_rlg,iyr,iset,zpdtim,
     &   oblat,oblon,obalt,asza,zenoff,azim,osds,
     &   opd,fovi,fovo,amal,ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,
     &   tins,pins,hins,tout,pout,hout,
     &   sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)

         write(*,'(a24,f8.3,f8.1,f9.3,f9.4,f10.6,2f8.1)')
     &   specname_cmc(:24),
     &   cleng(1),ctemp(1),cpres(1),pp(1,1),pp(1,1)/cpres(1),
     &   ifirst*graw,ilast*graw
         if(specname_cmc.ne.specname_rlg) then
             write(*,*) 'specname_cmc .ne. specname_rlg'
             write(*,*) specname_cmc//'  '//specname_rlg
             stop ' spectrum name mismatch'
         endif
         write(66,'(a24,f8.3,f8.1,f9.3,f9.4,f10.6,2f8.1)')
     &   specname_cmc(:24),
     &   cleng(1),ctemp(1),cpres(1),pp(1,1),pp(1,1)/cpres(1),
     &   ifirst*graw,ilast*graw

c Write a line of the .ray file.
         write(lunw_ray,format_ray) specname_cmc,zero,zero,zero,zero,
     &   zero,zero,(cleng(jcell)/1000.,jcell=1,ncell),
     &    (zero,jcell=ncell+1,mcell)

         write(lunw_mav,'(a)') 'Next Spectrum:'//specname_cmc
         write(lunw_mav,*) 2,4+nspeci_iso,ncell
         write(lunw_mav,*) header(:lh)
c
c  Write out the data lines of the .mav file.
         do jcell=1,ncell
            do jspeci=1,nspeci_iso
               vmr_mav(jspeci)=0.0
            end do
            tpp=0.0
            do igas=1,ngas(jcell)

               tpp=tpp+pp(igas,jcell)
               if(kiso_cmc(igas,jcell).eq.-1) then
                  do jspeci=1,nspeci_iso
c                     write(*,*) jcell,igas,jspeci,gasindex(jspeci),
c     &          kgas_cmc(igas,jcell)
                     if(gasindex(jspeci) .eq. kgas_cmc(igas,jcell)) then
                        vmr_mav(jspeci)=pp(igas,jcell)/cpres(jcell)
c                        write(*,*)jcell,igas,jspeci,pp(igas,jcell),cpres(jcell),
c     &                  vmr_mav(jspeci)
                  if(gasindex(jspeci).eq.49 .or. gasindex(jspeci).eq.71)
     &                  vmr_mav(jspeci)=vmr_mav(jspeci)/
     &                  fia(specindex(kgas_cmc(igas,jcell))+1)
                     endif
                  end do
               else
                  jspeci=specindex(kgas_cmc(igas,jcell))+
     &            kiso_cmc(igas,jcell)
                  vmr_mav(jspeci)=pp(igas,jcell)/cpres(jcell)/
     &            fia(jspeci)
               endif
c               write(*,*) specname_cmc(:12),jcell,igas,
c     &         kgas_cmc(igas,jcell),kiso_cmc(igas,jcell),
c     &         jspeci,pp(igas,jcell),vmr_mav(jspeci)
            end do   !  do igas=1,ngas(jcell)
c            write(*,*) (vmr_mav(specindex(6)+j),j=1,3)

c   Warn if sum of partial pressures exceeds total pressure.
            if(tpp.gt.cpres(jcell)) write(*,*)'Warning: SUM[PP] > Ptot',
     &      jcell,tpp,cpres(jcell),specname_cmc(:20)

            if(ctemp(jcell).lt.0.0) stop 'T < 0K'
            dens=.101325*cpres(jcell)/760./(ctemp(jcell)*1.38065E-23)
            write(lunw_mav,'(2f7.2,232(1pe11.3))')
     &      float(jcell),ctemp(jcell), cpres(jcell)/760.0,
     &      dens,(vmr_mav(jspeci),jspeci=1,nspeci_iso)
c            if(jcell.eq.1) write(*,*) 'slant col = ',
c     &       dens,vmr_mav(4),cleng(1),
c     &       dens*vmr_mav(4)*cleng(1)*100.0
         end do ! jcell=1,ncell

      end do  ! irec=1,999999
99    close(lunr_cmc)
      close(lunr_rlg)
      close(lunw_ray)
      close(lunw_mav)
      stop
      end
