c  summarize_linelist.f
c  Summarizes the content of a linelist in terms of the number of lines
c  the highest and lowest frequencies, strengths, and widths.
c
c  Produces three output files:
c  1) summary.rpt       A list of anomalies (e.g. -ve strengths, width=0)
c  2) summary_by_gas    A summary table on a per molecule basis
c  3) summary_by_speci  A summary table on a per speci (isotopolg) basis
c
      include "../ggg_int_params.f"

      integer mgas,ngas,lunr,lunw,luni,lun_rpt,
     & kgas,kiso,
     & jspeci,
     & nvmode
      parameter(lunr=14,lunw=15,luni=16,lun_rpt=17,
     & mgas=68)
      integer*4 count_gas(mgas),count_speci(mspeci),
     & isomax(mgas),molewt(mspeci),
     & tcount,istat,specindex(mgas+1),dgen(mvmode)
      character listnam*32,gasname*8,speci_id(mspeci)*24,
     & comment*120,gggdir*(mpath),dl*1
      real*8 freq,stren,r2,abhw,sbhw,eprime,td,pshift,
     & vmin_gas(mgas),vmax_gas(mgas),
     & vmin_speci(mspeci), vmax_speci(mspeci),
     & smin_gas(mgas),smax_gas(mgas),
     & smin_speci(mspeci),smax_speci(mspeci),
     & abhwmin_gas(mgas),abhwmax_gas(mgas),
     & abhwmin_speci(mspeci),abhwmax_speci(mspeci)
      real*4 fia,epsilon,
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
      data count_speci/mspeci*0/
      data vmin_speci/mspeci*0/
      data vmax_speci/mspeci*0/
      data smin_speci/mspeci*1/
      data smax_speci/mspeci*0/
      data abhwmin_speci/mspeci*1/
      data abhwmax_speci/mspeci*0/

      write(*,*)' summarize_linelist   version 1.0.3    30-Jul-2010'
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
      end do  ! jspeci=1,mspeci
      close(luni)
      specindex(kgas+1)=jspeci
      ngas=kgas
      nspeci=jspeci-1
      if(ngas.gt.mgas) STOP ' Increase parameter MGAS'

c============================================================
c  Read the selected linelist
      open(lun_rpt,file='anomaly.rpt', status='unknown')
      write(*,*)' Enter name of linelist (e.g. atm.101)'
      read(*,'(a)') listnam
      open(lunr,
     & file=gggdir(:lnbc(gggdir))//'/linelist/'//listnam,status='old')
      do iline=1,9999999 
      read(lunr,'(i2,i1,f12.6,2e10.3,2f5.4,f10.4,f4.2,f8.6,a)',end=33)
     & kgas,kiso,freq,stren,r2,abhw,sbhw,eprime,td,pshift,comment
      if(index(listnam,'hitran').gt.0) then
         call hitran_to_atmos_gas_numbering(kgas,kiso)
         if(kgas.le.0) then
            write(lun_rpt,*) kgas,kiso,freq,stren
            cycle
         endif
      endif
      kspeci= specindex(kgas)+kiso

c  Check that all linelist gases/species are in isotopomers.dat file
      if(kgas.le.0 .or. kgas.gt.ngas) then
         write(lun_rpt,*)' Warning: unknown gas:', kgas,kiso,freq
         cycle
      endif

      if(kspeci.gt.specindex(kgas+1)) then
         Write(lun_rpt,*)'Warning: unknown speci:',kgas,kiso,freq
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
      if(stren.lt.smin_speci(kspeci)) smin_speci(kspeci)=stren
      if(stren.gt.smax_speci(kspeci)) smax_speci(kspeci)=stren
      if(abhw.lt.abhwmin_speci(kspeci)) abhwmin_speci(kspeci)=abhw
      if(abhw.gt.abhwmax_speci(kspeci)) abhwmax_speci(kspeci)=abhw
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
     & //'  ABHWmin  ABHWmax'
      tcount=0
      do kgas=1,mgas
      write(lunw,'(2i4,i8,2f9.2,2e10.3,2f7.3)') kgas,isomax(kgas),
     & count_gas(kgas), vmin_gas(kgas),vmax_gas(kgas),
     & smin_gas(kgas),smax_gas(kgas),
     & abhwmin_gas(kgas),abhwmax_gas(kgas)
      tcount=tcount+count_gas(kgas)
      end do
      close(lunw)

c============================================================
c  Write out the summary by specie (isotopolog)
      open(lunw,file='summary_by_speci_'//listnam,status='unknown')
      write(lunw,*)'Gas Iso   Count    Vmin     Vmax     Smin      Smax'
     & //'  ABHWmin  ABHWmax'

      open(unit=luni,file=
     & gggdir(:lnbc(gggdir))//'/isotopologs/isotopologs.dat',
     & status='old')
      do jspeci=1,mspeci
         call read_isotop(luni,kgas,kiso,gasname,speci_id(jspeci),
     &   fia,delta,epsilon,molewt(jspeci),tdrpf(jspeci),
     &   vibfrq,dgen,nvmode,mvmode,istat)
         if(istat.ne.0) exit
         write(lunw,'(2i4,i8,2f9.2,2e10.3,2f7.3)') kgas,kiso,
     &   count_speci(jspeci), vmin_speci(jspeci),vmax_speci(jspeci),
     &   smin_speci(jspeci),smax_speci(jspeci),
     &   abhwmin_speci(jspeci),abhwmax_speci(jspeci)
      end do  ! jspeci=1,mspeci
      close(luni)
      close(lunw)
      write(*,'(a,i8,2f12.4)')' Classified Lines =  ',tcount
      stop
      end
