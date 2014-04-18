c  Program to generate a solar pseudo transmittance spectrum
c  from the solar linelist solar_merged_yyyymmdd.108
c  using the subroutine  solar_pts.f

      implicit none

      integer*4
     & i,ldot,lr,lnbc,
     & mcp,ncp, ! Number of spectral points to be computed
     & lunw,    ! LUN to write output file
     & lun_sll  ! Logical unit number for solar linelist
      parameter (lun_sll=14,lunw=15,mcp=2120000)

      character
     & outfile*80,
     & dl*1,
     & fmt*28,
     & solarll*80,  ! name of solar linelist
     & gggpath*180  ! gggpath

      real*4
     & spts(mcp)   !  Vector into which solar spectrum will be written

      real*8
     & frac,    ! Fraction of solar diameter viewed (0=disk center; 1=disk-integrated)
     & osvc,    ! Observer-Sun Velocity Component (m/s)
     & sol,     ! Speed of Light (m/s)
     & nus,     ! Starting frequency of interest (cm-1)
     & nue,     ! Ending frequency of interest (cm-1)
     & grid,    ! Frequency spacing out output spectrum (cm-1)
     & geff,    ! Effective grid point spacing (cm-1)
     & fzero    ! 

      parameter (sol=3.0E+08)  ! Speed of Light (m/s)

      write(*,*)' Enter starting and ending frequencies (cm-1):'
      read(*,*) nus, nue
      write(*,*)' Enter frequency point spacing (cm-1):'
      read(*,*) grid
      write(*,*)' Fraction of solar diameter observed '
      write(*,*)' (0=disc-center; 1=disk-integrated)'
      read(*,*) frac
      ncp=1+(nue-nus)/grid  ! number of spectral points
      if(ncp.gt.mcp) then
          write(*,*) 'ncp,mcp=',ncp,mcp
          stop ' Increase parameter MCP'
      endif

c      osvc=2500.0d0   ! Observer-Sun velocity component (m/s)
      osvc=0.0d0      ! Observer-Sun velocity component (m/s)

c      solarll='solar_merged_20121017.108'
      write(*,'(a)')
     & 'Enter the name of solar linelist (e.g. solar_merged.108): '
      read(*,'(a)')solarll

      call get_ggg_environment(gggpath,dl)
      lr = lnbc(gggpath)

      fzero=(nus-grid)*(1.d0+osvc/sol)
      geff=grid*(1.d0+osvc/sol)

      call solar_pts(lun_sll,gggpath(:lr)//dl//'linelist'//dl//solarll,
     & fzero,geff,frac,spts,ncp)

      if(nue .le. 9999.5) then
        fmt='(a1,i4.4,a1,i4.4,a1,i3.3,a4)'
      elseif(nus .ge. 9999.5) then
        fmt='(a1,i5.5,a1,i5.5,a1,i3.3,a4)'
      else
        fmt='(a1,i4.4,a1,i5.5,a1,i3.3,a4)'
      endif
   
      outfile=solarll
      ldot=index(outfile,'.')
      write(outfile(ldot:),fmt) '_',nint(nus),'_',nint(nue),'_',
     & nint(100*frac),'.out'
      open(lunw,file=outfile, status='unknown')
      do i=1,ncp
         write(lunw,'(f9.3,f9.5)') nus+(i-1)*grid,spts(i)
      end do
      close(lunw)
      stop
      end

c      include 'solar_pts.f'
c      include 'file_size_in_bytes.f'
c      include 'posnall.f'
c      include 'lnbc.f'
