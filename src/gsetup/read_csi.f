      subroutine read_csi(specname,trlg,lunr_csi,ncell,tins,
     & kcell,igas_in_cell,p_cell,t_cell,vmr_cell,cell_length)
c
c  Inputs:
c     specname            C**  Spectrun Name
c     trlg                R*8  ZPD Time
c     lunr_csi            I*4  LUN
c     ncell               I*4  Number of allowed cells
c     tins                R*8  Temperature inside instrument (C)
c
c  Outputs:
c     kcell               I*4  Current Cell #
c     igas_in_cell(ncell) I*4  Gas # in cells
c     p_cell(ncell)       I*4  Pressure in cells
c     t_cell(ncell)       I*4  Temp in cells
c     vmr_cell(ncell)     I*4  VMR in cell
c     cell_length(ncell)  I*4  Cell Length
c
c  Program to read the xx_cell_status_info.dat file and interpolate
c  conditions to trlg, the ZPD time of the spectrum in the runlog.
c  Code figures out whether there was a cell in use, or not. And
c  figures out when to read another record from cell_status_info.dat.
c
c  In the event that the runlog contains spectra from multiple sites,
c  the relevant xx_cell_status_info.dat files are closed and opened.
c
c  Can handle up to 2 different cell types per site.
c
c  Assumes that:
c  1) Runlog is chronologically ordered
c  2) cell status files are chronologically ordered
c  3) Only one cell is present at a time
c
c  The cell lengths computed here are effective cell lengths,
c  not the true cell lengths. This is because Frank Hase measured
c  effective cell pressures defined by
c   Peff.SBHW_hcl = PP_hcl.SBHW_hcl + PP_air.ABHW_hcl + PP_h2o.WBHW_hcl
c  so that
c   Peff = PP_hcl + PP_air.ABHW_hcl/SBHW_hcl + PP_h2o.WBHW_hcl/SBHW_hcl
c  Peff is the pressure of a pure HCl sample that would produce the
c  same line broadening as was measured. Since the gas cell contains
c  contaminants (e.g. air, H2O) which also broaden the HCl lines,
c  Peff is higher than the true PP_hcl. The use of effective pressure
c  requires the assumption that the HCl vmr inside the cell is 1.0,
c  otherwise GFIT will compute additional broadening from air. And
c  given that the VMR=1, the only way to match the slant column is
c  to compute an effective cell length such that
c      Cell_Length_eff = Column.Boltzman.T_hcl / P_eff.
c  where T is the temperature (K) of the HCl. So this explains why
c  there is no cell length information in the xx_cell_status_info.dat
c  files. Even if it was there, we wouldn't be able to use it.
c
c
c Pseudo-code
c   
c Open relevant csi-file first time or when site changes
c     sp=site_prefix(specname)
c     if(sp.ne.spwas) then
c        close(lunr_csi)
c        site_prefix=specname(1:2)
c        open(lunr_csi,file=csi_file,status='old')
c        t2=-1.E+12  !  Long time ago
c        t2was=t2    !  Long time ago
c        cell_flag=0
c        kcell=1
c        do jcell=1,ncell
c           igas_in_cell(jcell)=0
c           p_cell(jcell)=0.0
c           t_cell(jcell)=tins
c           vmr_cell(jcell)=0.0
c           cell_length(jcell)=0.0
c        end do
c        spwas=sp
c     endif
c
c Read xx_cell_status_info.dat until it brackets trlg
c     do while (trlg.gt.t2 .and. cell_flag.ge.0)
c        t2was=t2
c        read(lunr_csi,csi_fmt,iostat=cell_flag) t1,t2,...
c        if(t1.lt.t2was) stop 'cell file not chronological?'
c     end do  ! while (trlg.gt.t2)
c
c Find KCELL that corresponds to KGAS.  If KGAS is new, create index.
c     kcellwas=kcell
c     do kcell=1,ncell
c        if(igas_in_cell(kcell).eq.0) then ! KCELL not yet indexed to KGAS
c           igas_in_cell(kcell)=kgas  ! New gas
c           exit
c        elseif(kgas.eq.igas_in_cell(kcell) ) then ! already indexed gas
c           exit
c        endif
c     end do ! kcell=1,ncell

c See if trlg corresponds to a time when cell was present, or not
c If cell present compute cell conditions. If not, set cell_length=0
c     if(trlg.gt.t2was) then  ! Cell removed
c        cell_length(kcellwas)=0.0
c     endif
c     if(trlg.gt.t2) then  ! No Cell. Past EOF of cell_status_info.dat
c        cell_length(kcell)=0.0
c     elseif(trlg.gt.t1) then  ! Cell
c        igas_in_cell(kcell)=kgas
c        vmr_cell(kcell)=1.0
c        p_cell(kcell)=
c        cell_length(kcell)=
c     elseif(trlg.lt.t2was) then   ! Something is wrong
c        stop 'Error (runlog not chronological?)'
c     endif

      implicit none
      integer*4 ncell,kcell,jcell,kcellwas,
     & lunr_csi,nlhead,ncol,i,kgas,
     & lnbc,
     & jd1,jd2,
     & cell_flag,
     & y1,m1,d1,h1,n1,s1,
     & y2,m2,d2,h2,n2,s2

      integer*4 igas_in_cell(ncell)

      real*4 boltzman,column,
     & p_cell(ncell),t_cell(ncell),vmr_cell(ncell),cell_length(ncell)

      real*8 column1,column2,peff1,peff2,tins
      real*8 t1,t2,t2was,trlg
      parameter (boltzman=1.38065E-23)

      character
     & csi_fmt*128,root*128,dl*1,csi_file*256,
     & specname*(*),sp*2,spwas*2,site_prefix*2,string*128

      save 
     & column1,column2,peff1,peff2,spwas,cell_flag,
     & t1,t2,t2was,csi_fmt,kgas,kcellwas

      data spwas/'  '/

      sp=site_prefix(specname)
      if(sp.ne.spwas) then  ! New Site
          close(lunr_csi)
         call get_ggg_environment(root, dl)
         csi_file=root(:lnbc(root))//'cell_status_info'//dl//
     &   sp//'_cell_status_info.dat'
         write(*,*)'Opening CSI file: '//csi_file(:60)
         open(lunr_csi,file=csi_file,status='old')
         read(lunr_csi,*) nlhead,ncol
         do i=2,nlhead
            read(lunr_csi,'(a)') string
            if(string(1:7).eq.'format=') csi_fmt=string(8:)
         end do
c         write(*,*)'csi_fmt = '//csi_fmt
         t2=-1.E+12  !  Long time ago
         t2was=t2    !  Long time ago
         cell_flag=0
         kcell=1
         do jcell=1,ncell
            igas_in_cell(jcell)=0
            p_cell(jcell)=0.0
            t_cell(jcell)=273.15+sngl(tins)
            vmr_cell(jcell)=0.0
            cell_length(jcell)=0.0
         end do
         spwas=sp
      endif

c  Read xx_cell_status_info.dat until it brackets trlg
      do while (trlg.gt.t2 .and. cell_flag.ge.0) 
         t2was=t2
         read(lunr_csi,csi_fmt,iostat=cell_flag)
     &   y1,m1,d1,h1,n1,s1,y2,m2,d2,h2,n2,s2,kgas,column1,column2,
     &   peff1,peff2
c         write(82,'(a,i5,6i3,4e12.4,i5)')'Read cell info: ',
c     &   y1,m1,d1,h1,n1,s1,kgas,column1,column2,peff1,peff2,cell_flag
         if(cell_flag.lt.0) exit  ! EOF

         call julian(y1,m1,d1,jd1)
         t1=jd1+(h1+(n1+dfloat(s1)/60)/60)/24
         call julian(y2,m2,d2,jd2)
         t2=jd2+(h2+(n2+dfloat(s2)/60)/60)/24
c          write(*,*) 't2 = ',y2,m2,d2,h2,n2,s2,t2
         if(t1.lt.t2was) stop 'cell file not chronological?'
      end do  ! while (trlg.gt.t2)
c
c  Find KCELL that corresponds to KGAS.  If KGAS is new, create index.
      kcellwas=kcell
      do kcell=1,ncell
         if(igas_in_cell(kcell).eq.0) then ! KCELL not yet indexed to KGAS
            igas_in_cell(kcell)=kgas  ! New gas
            exit
         elseif(kgas.eq.igas_in_cell(kcell) ) then ! already indexed gas
            exit
         endif
      end do ! kcell=1,ncell
c
      t_cell(kcell)=sngl(tins)+273.15

c  See if trlg corresponds to a time when cell was present, or not
      if(trlg.gt.t2was) then  ! Cell removed
         cell_length(kcellwas)=0.0
c         write(*,*)'kcell2',kcellwas,kcell,kgas
      endif

      if(trlg.gt.t2) then  ! No Cell. Past EOF of cell_status_info.dat
         cell_length(kcell)=0.0
c         write(*,*)'kcell0',kcellwas,kcell
      elseif(trlg.gt.t1) then  ! Cell
         igas_in_cell(kcell)=kgas
         vmr_cell(kcell)=1.0
         p_cell(kcell)=sngl(peff1+(trlg-t1)*(peff2-peff1)/(t2-t1))
         column=sngl(column1+(trlg-t1)*(column2-column1)/(t2-t1))
         cell_length(kcell)=column*boltzman*t_cell(kcell)/
     &   (100*p_cell(kcell))
         cell_length(kcell)=0.001*cell_length(kcell)  ! convert m to km
c         write(*,*)'kcell1',kcellwas,kcell,kgas
      elseif(trlg.lt.t2was) then       ! Something is wrong
         write(*,'(a20,4f15.5)') specname,trlg,t2was,t1,t2
         stop 'Error (runlog not chronological?)'
      endif
      return
      end

      function site_prefix(specname)
c  Derives a 2-character prefix from each spectrum name.
c  Necessary to handle non-TCCON spectrum names.
c  For example, prevents an ATMOS spectrum beginning "pat"
c  from being mistaken for a Park Falls TCCON spectrum.

      character site_prefix*2,specname*(*)

      if(specname(1:3).eq.'sao') then
         site_prefix='sa'           ! SAO
      elseif(specname(1:5).eq.'gosat') then
         site_prefix='go'           ! GOSAT
      elseif(specname(1:3).eq.'fsu') then
         site_prefix='fs'           ! 
      elseif(specname(1:3).eq.'irr') then
         site_prefix='ir'           ! 
      elseif(specname(1:3).eq.'ace') then
         site_prefix='ac'           ! ACE
      elseif(specname(1:2).eq.'ss') then
         site_prefix='ac'           ! ACE
      elseif(specname(1:2).eq.'du') then
         site_prefix='du'           ! Denver University
      elseif(specname(1:3).eq.'psl') then
         site_prefix='at'           ! ATMOS SL3
      elseif(specname(1:3).eq.'pat') then
         site_prefix='at'           ! ATMOS ATLAS
      elseif(specname(2:3).eq.'hg') then
         site_prefix='m4'           ! MkIV HgCdTe
      elseif(specname(2:3).eq.'in') then
         site_prefix='m4'           ! MkIV InSb
      elseif(abs(ichar(specname(2:2))-52.5).lt.5. ! digit
     & .and. specname(1:1).eq.'r') then
         site_prefix='at'           ! ATMOS 
      elseif(abs(ichar(specname(1:1))-52.5).lt.5. ! digit
     & .and. specname(7:7).eq.'R') then
         site_prefix='kp'           ! Kitt Peak
      elseif(specname(3:3).eq.'2') then
         site_prefix=specname(1:2)  ! TCCON
      else
         write(*,*)'site_prefix: '//specname
         stop ' unrecognised spectrum name'
      endif
      return
      end
      
