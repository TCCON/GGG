c  summarize_linelist.f
c  Summarizes the content of a linelist in terms of the number of lines
c  the highest and lowest frequencies, strengths, and widths.
c
c  Produces five output files:
c  1) summary.rpt       A list of anomalies (e.g. -ve strengths, width=0, duplicates)
c  2) summary_by_gas    A summary table on a per molecule basis
c  3) summary_by_speci  A summary table on a per speci (isotopolg) basis
c  And in the event of a mis-sorted linelist
c  4) ll_mis-sorted.161
c  5) ll_sorted.161
c
      implicit none
      include "../gfit/ggg_int_params.f"

      integer ngas,lunr_ll,lunw,lunr_iso,lunw_rpt,
     & lunw_cs,lunw_ms,ciso2kiso,ncs,nms,nlines,
     & i,igas,kgas,kiso,iiso,ll,idum,ls,ncpl,
     & jspeci,iline,kspeci,nspeci,lnblnk,
     & nvmode,nlhead,ncol,icode,ianom,idupl

      parameter(lunr_ll=14,lunw=15,lunr_iso=16,lunw_rpt=17,
     & lunw_cs=20, lunw_ms=21)

      integer*4 count_gas(mgas),count_speci(mspeci),
     & isomin(mgas),isomax(mgas),molewt(mspeci),
     & tcount,istat,specindex(mgas+1),dgen(mvmode)

      integer*8 fsib,file_size_in_bytes

      character listnam*80,gasname*8,speci_id(mspeci)*20,
     & iso_fmt*80,ciso*1,cpl*3,
     & llstring*240,llswas*240,comment*93,gggdir*(mpath),dl*1,version*56

      real*8 fwas,freq,stren,r2,abhw,sbhw,eprime,td,pshift,
     & vmin_gas(mgas),vmax_gas(mgas),
     & vmin_speci(mspeci), vmax_speci(mspeci),
     & smin_gas(mgas),smax_gas(mgas),
     & smin_speci(mspeci),smax_speci(mspeci),
     & abhwmin_gas(mgas),abhwmax_gas(mgas),
     & abhwmin_speci(mspeci),abhwmax_speci(mspeci),
     & sbhwmin_gas(mgas),sbhwmax_gas(mgas),
     & sbhwmin_speci(mspeci),sbhwmax_speci(mspeci),
     & pshmin_gas(mgas),pshmax_gas(mgas),
     & pshmin_speci(mspeci),pshmax_speci(mspeci)
      real*4 fia,lnfrd,delta,
     & tdrpf(mspeci),ewvb(mspeci),atc(mspeci),vibfrq(mvmode)
c============================================================
      data isomin/mgas*99/
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
      data pshmin_gas/mgas*1/
      data pshmax_gas/mgas*-1/
      data count_speci/mspeci*0/
      data vmin_speci/mspeci*0/
      data vmax_speci/mspeci*0/
      data smin_speci/mspeci*1/
      data smax_speci/mspeci*0/
      data abhwmin_speci/mspeci*1/
      data abhwmax_speci/mspeci*0/
      data sbhwmin_speci/mspeci*1/
      data sbhwmax_speci/mspeci*0/
      data pshmin_speci/mspeci*1/
      data pshmax_speci/mspeci*-1/

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

      version=' summarize_linelist       version 1.22       2020-06-20 '
      write(*,*) version
      write(*,*)
      call get_ggg_environment(gggdir, dl)
      ianom=0
      idupl=0
c============================================================
c  Read isotopologs.dat file
      open(unit=lunr_iso,file=
     & gggdir(:lnblnk(gggdir))//'isotopologs'//dl//'isotopologs.dat',
     & status='old')
      read(lunr_iso,*) nlhead,ncol
      read(lunr_iso,'(7x,a)') iso_fmt
      do i=3,nlhead
         read(lunr_iso,*)
      end do
      do jspeci=1,mspeci
         call read_isotopolog(lunr_iso,iso_fmt,igas,iiso,gasname,
     &   speci_id(jspeci),icode,fia,delta,lnfrd,molewt(jspeci),
     &   ewvb(jspeci),atc(jspeci),tdrpf(jspeci),vibfrq,dgen,nvmode,
     &   mvmode,istat)
         if(istat.ne.0) exit
c         write(*,*) igas, iiso, gasname, speci_id(jspeci)
         if(igas.lt.0) stop 'KGAS<0'
         if(iiso.lt.0) stop 'KISO<0'
         specindex(igas)=jspeci-iiso
c         write(*,*)jspeci,igas,iiso,specindex(igas)
      end do  ! jspeci=1,mspeci
      close(lunr_iso)
      specindex(igas+1)=jspeci
      ngas=igas
      nspeci=jspeci-1
      if(ngas.gt.mgas) STOP ' Increase parameter MGAS'
c============================================================
c  Read the selected linelist
      if (iargc() == 0) then
         write(*,*)' Enter name of linelist (e.g. atm.161)'
         read(*,'(a)') listnam
      elseif (iargc() == 1) then
         call getarg(1, listnam)
      else
         stop 'Usage: $gggpath/bin/summarize_linelist '//
     & 'linelistname (e.g. atm.161)'
      endif
      ll=lnblnk(listnam)
      cpl=listnam(ll-2:ll)
      read(cpl,*) ncpl
      if(ncpl.eq.108) stop 'unrecognized linelist format'
      fsib=file_size_in_bytes(lunr_ll,listnam)
      nlines=int(fsib/ncpl,kind(ncpl))
      if ( nlines*ncpl .ne. fsib ) then
         write(*,*)'Warning: Linelist size not exactly divisible by '//
     &   'record length'
         write(*,*)'Linelist size = ',fsib
         write(*,*)'Record Length = ',ncpl
      endif

      open(lunr_ll,
c     & file=gggdir(:lnblnk(gggdir))//'/linelist/'//listnam,status='old')
     & file=listnam,status='old')
      open(lunw_rpt,file='summarize_linelist_anomalies',
     & status='replace')
      ncs=0
      open(lunw_cs,file='sorted_lines.'//cpl,
     & status='replace')
      nms=0
      open(lunw_ms,file='mis-sorted_lines.'//cpl,
     & status='replace')
      fwas=0.0d0
      do iline=1,9999999 
         read(lunr_ll,'(a)',end=33) llstring
         ls=lnblnk(llstring)
         if(ichar(llstring(ncpl+1:ncpl+1)).ne.32) 
     &   write(*,*) iline,ichar(llstring(ncpl+1:ncpl+1))
         if(ls.le.0) then
            write(lunw_rpt,*)iline,' blank record'
            if(ichar(llstring(ncpl:ncpl)).ne.32) write(*,*) iline
            ianom=ianom+1
         endif
         read(llstring,'(i2,a1,f12.6,2e10.3,2f5.4,f10.4,f4.2,f8.6,a)')
     &    kgas,ciso,freq,stren,r2,abhw,sbhw,eprime,td,pshift,comment
         kiso=ciso2kiso(kgas,ciso)
c         write(*,*) kgas,kiso,freq,stren

         if(llstring.eq.llswas) then
            write(lunw_rpt,'(a,i8,2i3,f11.4,e10.3)')
     &      'consecutive duplicated lines',iline,kgas,kiso,freq,stren
            ianom=ianom+1
            idupl=idupl+1
         endif
             
         if(freq.gt.fwas .or.
     &      (abs(freq-fwas).le.0.0d0 .and. lle(llswas,llstring))) then
            write(lunw_cs,'(a)') llstring(:ncpl-1)  ! Correctly sorted
            ncs=ncs+1
         else
            write(lunw_rpt,'(a,i7,2i3,2f11.4,f8.4)')
     &      'freq<=fwas:',iline,kgas,kiso,fwas,freq,freq-fwas
            write(lunw_ms,'(a)') llstring(:ncpl-1)  ! Mis-sorted
            nms=nms+1
         endif
         fwas=freq

c         write(*,*) iline,ls,kgas,kiso,freq,stren
         if(index(listnam,'hitran').gt.0) then
            call hitran_to_atmos_gas_numbering(kgas,kiso)
            if(kgas.le.0) cycle
         endif
         if (kgas.eq.0) write(*,*)'Error: kgas=0  linelist line #',iline
         kspeci= specindex(kgas)+kiso

c  Check that all linelist gases/species are in isotopomers.dat file
         if(kgas.le.0 .or. kgas.gt.ngas) then
            write(lunw_rpt,*)'Warning:unknown gas:',iline,kgas,kiso,freq
            ianom=ianom+1
            cycle
         endif

         if(kspeci.gt.specindex(kgas+1)) then
            write(lunw_rpt,'(a,3i8,f12.6,2i8)')'Warning:unknown speci:',
     &      iline,kgas,kiso,freq,kspeci,specindex(kgas+1)
            ianom=ianom+1
            cycle
         endif

         if(stren.le.0 .or. abhw.le.0.0 .or. sbhw.le.0.0) then
            if(kgas.ne.54 .and. kgas.ne.72) then ! skip Mars & Venus dust
               write(lunw_rpt,'(i2,a1,f12.6,e10.3,2(1x,f5.4))')
     &         kgas,ciso,freq,stren,abhw,sbhw
               ianom=ianom+1
            endif
         endif
         if(stren.lt.smin_gas(kgas)) smin_gas(kgas)=stren
         if(stren.gt.smax_gas(kgas)) smax_gas(kgas)=stren
         if(abhw.lt.abhwmin_gas(kgas)) abhwmin_gas(kgas)=abhw
         if(abhw.gt.abhwmax_gas(kgas)) abhwmax_gas(kgas)=abhw
         if(sbhw.lt.sbhwmin_gas(kgas)) sbhwmin_gas(kgas)=sbhw
         if(sbhw.gt.sbhwmax_gas(kgas)) sbhwmax_gas(kgas)=sbhw
         if(pshift.lt.pshmin_gas(kgas)) pshmin_gas(kgas)=pshift
         if(pshift.gt.pshmax_gas(kgas)) pshmax_gas(kgas)=pshift

         if(stren.lt.smin_speci(kspeci)) smin_speci(kspeci)=stren
         if(stren.gt.smax_speci(kspeci)) smax_speci(kspeci)=stren
         if(abhw.lt.abhwmin_speci(kspeci)) abhwmin_speci(kspeci)=abhw
         if(abhw.gt.abhwmax_speci(kspeci)) abhwmax_speci(kspeci)=abhw
         if(sbhw.lt.sbhwmin_speci(kspeci)) sbhwmin_speci(kspeci)=sbhw
         if(sbhw.gt.sbhwmax_speci(kspeci)) sbhwmax_speci(kspeci)=sbhw
         if(pshift.lt.pshmin_speci(kspeci)) pshmin_speci(kspeci)=pshift
         if(pshift.gt.pshmax_speci(kspeci)) pshmax_speci(kspeci)=pshift

         if(count_gas(kgas).eq.0) vmin_gas(kgas)=freq
         if(count_speci(kspeci).eq.0) vmin_speci(kspeci)=freq
         count_gas(kgas)=count_gas(kgas)+1
         count_speci(kspeci)=count_speci(kspeci)+1
         if(kiso.lt.isomin(kgas)) isomin(kgas)=kiso
         if(kiso.gt.isomax(kgas)) isomax(kgas)=kiso
         vmax_gas(kgas)=freq
         vmax_speci(kspeci)=freq
         if(mod(iline,100000).eq.0) write(*,*) iline, ' lines'
         llswas=llstring
      end do
33    close(lunr_ll)
      close(lunw_rpt)
      close(lunw_cs)
      close(lunw_ms)
      write(*,*) 'Total number of Lines =  ',iline-1
      write(*,*) 'Correctly-Sorted Lines = ',ncs
      write(*,*) 'Mis-Sorted Lines =       ',nms
      write(*,*) 'Duplicated Lines =       ',idupl

c============================================================
c  Write out the summary by gas
      open(lunw,file='summary_by_gas_'//listnam(:ll-4),status='unknown')
      write(lunw,*)3,14
      write(lunw,*)version
      write(lunw,*)'Gas Imin Imax  Count    Vmin     Vmax    Smin     
     &Smax   ABHWmin ABHWmax SBHWmin SBHWmax PShmin PShmax'
      tcount=0
      do kgas=1,mgas
         write(lunw,'(3i4,i8,2f9.2,2(1pe10.3),4(0pf8.3),2f8.4)') kgas,
     &    isomin(kgas),isomax(kgas),
     &    count_gas(kgas), vmin_gas(kgas),vmax_gas(kgas),
     &    smin_gas(kgas),smax_gas(kgas),
     &    abhwmin_gas(kgas),abhwmax_gas(kgas),
     &    sbhwmin_gas(kgas),sbhwmax_gas(kgas),
     &    pshmin_gas(kgas),pshmax_gas(kgas)
         tcount=tcount+count_gas(kgas)
      end do
      close(lunw)
      write(*,*) 'Classifiable Lines =     ',tcount
      write(*,*) 'Unclassifiable Lines =   ',iline-1-tcount
c============================================================
c  Write out the summary by specie (isotopolog)
      open(lunw,file='summary_by_speci_'//listnam(:ll-4),
     & status='unknown')
      write(lunw,*)3,13
      write(lunw,*)version
      write(lunw,*)'Gas Iso   Count    Vmin     Vmax     Smin      Smax'
     & //'    ABHWmin ABHWmax SBHWmin SBHWmax PShmin PShmax' 

      open(unit=lunr_iso,file=
     & gggdir(:lnblnk(gggdir))//'isotopologs'//dl//'isotopologs.dat',
     & status='old')
      read(lunr_iso,*) nlhead,ncol
      read(lunr_iso,'(7x,a)') iso_fmt
      do i=3,nlhead
         read(lunr_iso,*)
      end do
      do jspeci=1,mspeci
         call read_isotopolog(lunr_iso,iso_fmt,igas,iiso,gasname,
     &   speci_id(jspeci),icode,fia,delta,lnfrd,molewt(jspeci),
     &   ewvb(jspeci),atc(jspeci),tdrpf(jspeci),vibfrq,dgen,nvmode,
     &   mvmode,istat)
         if(istat.ne.0) exit
         write(lunw,'(2i4,i8,2f9.2,2(1pe10.3),4(0pf8.3),2f8.4)')
     &   igas,iiso,
     &   count_speci(jspeci), vmin_speci(jspeci),vmax_speci(jspeci),
     &   smin_speci(jspeci),smax_speci(jspeci),
     &   abhwmin_speci(jspeci),abhwmax_speci(jspeci),
     &   sbhwmin_speci(jspeci),sbhwmax_speci(jspeci),
     &   pshmin_speci(jspeci),pshmax_speci(jspeci)
      end do  ! jspeci=1,mspeci
      close(lunr_iso)
      close(lunw)
      write(*,'(a,i8)')' Anomalies =  ',ianom
      stop
      end
