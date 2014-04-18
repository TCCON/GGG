c  summarize_linelist.f
c  Summarizes the content of a linelist in terms of the number of lines
c  the highest and lowest frequencies, strengths, and widths.
c
c  Produces three output files:
c  1) summary.rpt       A list of anomalies (e.g. -ve strengths, width=0)
c  2) summary_by_gas    A summary table on a per molecule basis
c  3) summary_by_speci  A summary table on a per speci (isotopolg) basis
c
      implicit none
      include "../ggg_int_params.f"

      integer ngas,lunr,lunw,luni,lun_rpt,
     & kgas,kiso,lnbc,ll,
     & jspeci,iline,kspeci,nspeci,
     & nvmode
      parameter(lunr=14,lunw=15,luni=16,lun_rpt=17)
      integer*4 count_gas(mgas),count_speci(mspeci),
     & isomax(mgas),molewt(mspeci),
     & tcount,istat,specindex(mgas+1),dgen(mvmode)
      character listnam*80,gasname*8,speci_id(mspeci)*24,
     & string*240,comment*120,gggdir*(mpath),dl*1,version*48
      real*8 freq,stren,r2,abhw,sbhw,eprime,td,pshift,
     & vmin_gas(mgas),vmax_gas(mgas),
     & vmin_speci(mspeci), vmax_speci(mspeci),
     & smin_gas(mgas),smax_gas(mgas),
     & smin_speci(mspeci),smax_speci(mspeci),
     & abhwmin_gas(mgas),abhwmax_gas(mgas),
     & abhwmin_speci(mspeci),abhwmax_speci(mspeci),
     & sbhwmin_gas(mgas),sbhwmax_gas(mgas),
     & sbhwmin_speci(mspeci),sbhwmax_speci(mspeci)
      real*4 fia,epsilon,delta,
     & tdrpf(mspeci),vibfrq(mvmode)
c============================================================
      data isomax/mgas*0/
      data count_gas/mgas*0/
      data vmin_gas/mgas*0/
      data vmax_gas/mgas*0/
      data smin_gas/mgas*1/
      data smax_gas/mgas*0/
      data abhwmin_gas/mgas*1/
      data abhwmax_gas/mgas*0/
      data sbhwmin_gas/mgas*1/
      data sbhwmax_gas/mgas*0/
      data count_speci/mspeci*0/
      data vmin_speci/mspeci*0/
      data vmax_speci/mspeci*0/
      data smin_speci/mspeci*1/
      data smax_speci/mspeci*0/
      data abhwmin_speci/mspeci*1/
      data abhwmax_speci/mspeci*0/
      data sbhwmin_speci/mspeci*1/
      data sbhwmax_speci/mspeci*0/

      version=' summarize_linelist   version 1.08   2013-05-22 '
      write(*,*) version
      write(*,*)
      call get_ggg_environment(gggdir, dl)

c============================================================
c  Read isotopologs.dat file
      open(unit=luni,file=
     & gggdir(:lnbc(gggdir))//'isotopologs'//dl//'isotopologs.dat',
     & status='old')
      do jspeci=1,mspeci
         call read_isotop(luni,kgas,kiso,gasname,speci_id(jspeci),
     &   fia,delta,epsilon,molewt(jspeci),tdrpf(jspeci),
     &   vibfrq,dgen,nvmode,mvmode,istat)
         if(istat.ne.0) exit
c         write(*,*) kgas, kiso, gasname, speci_id(jspeci)
         if(kgas.lt.0) stop 'KGAS<0'
         if(kiso.lt.0) stop 'KISO<0'
         specindex(kgas)=jspeci-kiso
         write(*,*)jspeci,kgas,kiso,specindex(kgas)
      end do  ! jspeci=1,mspeci
      close(luni)
      specindex(kgas+1)=jspeci
      ngas=kgas
      nspeci=jspeci-1
      if(ngas.gt.mgas) STOP ' Increase parameter MGAS'

c============================================================
c  Read the selected linelist
      if (iargc() == 0) then
         write(*,*)' Enter name of linelist (e.g. atm.101)'
         read(*,'(a)') listnam
      elseif (iargc() == 1) then
         call getarg(1, listnam)
      else
         stop 'Usage: $gggpath/bin/summarize_linelist '//
     & 'linelistname (e.g. atm.101)'
      endif
      ll=lnbc(listnam)
      if(listnam(ll-3:ll).eq.'.108') then
        stop 'unrecognized linelist format'
      endif
      open(lunr,
c     & file=gggdir(:lnbc(gggdir))//'/linelist/'//listnam,status='old')
     & file=listnam,status='old')
      open(lun_rpt,file='anomaly_'//listnam,
     & status='unknown')
      do iline=1,9999999 
      read(lunr,'(a)',end=33) string
      if(lnbc(string).le.0) write(lun_rpt,*)iline,' blank record'
      if(ichar(string(161:161)).ne.32) write(*,*) iline
      if(ichar(string(162:162)).ne.32) write(*,*) iline
      read(string,'(i2,i1,f12.6,2e10.3,2f5.4,f10.4,f4.2,f8.6,a)')
     & kgas,kiso,freq,stren,r2,abhw,sbhw,eprime,td,pshift,comment
      if(kgas.eq.2 .and. kiso.eq.0) kiso=10  ! HITRAN 2012 kludge
c      write(*,*) kgas,kiso,freq,stren
      if(index(listnam,'hitran').gt.0) then
         call hitran_to_atmos_gas_numbering(kgas,kiso)
         if(kgas.le.0) then
c            write(lun_rpt,*) iline,kgas,kiso,freq,stren
            cycle
         endif
      endif
      kspeci= specindex(kgas)+kiso

c  Check that all linelist gases/species are in isotopomers.dat file
      if(kgas.le.0 .or. kgas.gt.ngas) then
         write(lun_rpt,*)' Warning: unknown gas:', iline,kgas,kiso,freq
         cycle
      endif

      if(kspeci.gt.specindex(kgas+1)) then
         write(lun_rpt,'(a,3i8,f12.6,2i8)')'Warning: unknown speci:',
     &   iline,kgas,kiso,freq,kspeci,specindex(kgas+1)
         cycle
      endif

      if(stren.le.0 .or. abhw.le.0.0)
     &  write(lun_rpt,'(i2,i1,f12.6,e10.3,2(1x,f5.4))')
     & kgas,kiso,freq,stren,abhw,sbhw
c      if(abhw.le.0.0 .and. kgas.eq.30)
c     & write(25,'(i2,i1,f12.6,2e10.3,2f5.4,f10.4,f4.2,f8.6,a)')
c     & kgas,kiso,freq,stren,r2,abhw,sbhw,eprime,td,pshift,comment
      if(stren.lt.smin_gas(kgas)) smin_gas(kgas)=stren
      if(stren.gt.smax_gas(kgas)) smax_gas(kgas)=stren
      if(abhw.lt.abhwmin_gas(kgas)) abhwmin_gas(kgas)=abhw
      if(abhw.gt.abhwmax_gas(kgas)) abhwmax_gas(kgas)=abhw
      if(sbhw.lt.sbhwmin_gas(kgas)) sbhwmin_gas(kgas)=sbhw
      if(sbhw.gt.sbhwmax_gas(kgas)) sbhwmax_gas(kgas)=sbhw
      if(stren.lt.smin_speci(kspeci)) smin_speci(kspeci)=stren
      if(stren.gt.smax_speci(kspeci)) smax_speci(kspeci)=stren
      if(abhw.lt.abhwmin_speci(kspeci)) abhwmin_speci(kspeci)=abhw
      if(abhw.gt.abhwmax_speci(kspeci)) abhwmax_speci(kspeci)=abhw
      if(sbhw.lt.sbhwmin_speci(kspeci)) sbhwmin_speci(kspeci)=sbhw
      if(sbhw.gt.sbhwmax_speci(kspeci)) sbhwmax_speci(kspeci)=sbhw
      if(count_gas(kgas).eq.0) vmin_gas(kgas)=freq
      if(count_speci(kspeci).eq.0) vmin_speci(kspeci)=freq
      count_gas(kgas)=count_gas(kgas)+1
      count_speci(kspeci)=count_speci(kspeci)+1
      if(kiso.gt.isomax(kgas)) isomax(kgas)=kiso
      vmax_gas(kgas)=freq
      vmax_speci(kspeci)=freq
      if(mod(iline,100000).eq.0) write(*,*) iline, ' lines'
      end do
33    close(lunr)
      close(lun_rpt)
      write(*,*) 'Total lines = ',iline-1

c============================================================
c  Write out the summary by gas
      open(lunw,file='summary_by_gas_'//listnam,status='unknown')
      write(lunw,*)'Gas Imax  Count    Vmin     Vmax     Smin      Smax'
     & //'    ABHWmin ABHWmax SBHWmin SBHWmax'
      tcount=0
      do kgas=1,mgas
      write(lunw,'(2i4,i8,2f9.2,2e10.3,4f8.3)') kgas,isomax(kgas),
     & count_gas(kgas), vmin_gas(kgas),vmax_gas(kgas),
     & smin_gas(kgas),smax_gas(kgas),
     & abhwmin_gas(kgas),abhwmax_gas(kgas),
     & sbhwmin_gas(kgas),sbhwmax_gas(kgas)
      tcount=tcount+count_gas(kgas)
      end do
      close(lunw)

c============================================================
c  Write out the summary by specie (isotopolog)
      open(lunw,file='summary_by_speci_'//listnam,status='unknown')
      write(lunw,*)'Gas Iso   Count    Vmin     Vmax     Smin      Smax'
     & //'    ABHWmin ABHWmax SBHWmin SBHWmax' 

      open(unit=luni,file=
     & gggdir(:lnbc(gggdir))//'/isotopologs/isotopologs.dat',
     & status='old')
      do jspeci=1,mspeci
         call read_isotop(luni,kgas,kiso,gasname,speci_id(jspeci),
     &   fia,delta,epsilon,molewt(jspeci),tdrpf(jspeci),
     &   vibfrq,dgen,nvmode,mvmode,istat)
         if(istat.ne.0) exit
         write(lunw,'(2i4,i8,2f9.2,2e10.3,4f8.3)') kgas,kiso,
     &   count_speci(jspeci), vmin_speci(jspeci),vmax_speci(jspeci),
     &   smin_speci(jspeci),smax_speci(jspeci),
     &   abhwmin_speci(jspeci),abhwmax_speci(jspeci),
     &   sbhwmin_speci(jspeci),sbhwmax_speci(jspeci)
      end do  ! jspeci=1,mspeci
      close(luni)
      close(lunw)
      write(*,'(a,i8,2f12.4)')' Classified Lines =  ',tcount
      stop
      end
