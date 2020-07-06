c  Program to generate a solar pseudo transmittance spectrum
c  from the solar linelist solar_merged_yyyymmdd.108
c  using the subroutine  solar_pts.f

      implicit none

      integer*4
     & i,ldot,lr,lnbc,
     & mcp,ncp, ! Number of spectral points to be computed
     & ifsp,    ! Index of First Spectral Point
     & lunw,    ! LUN to write output file
     & lun_sll  ! Logical unit number for solar linelist
      parameter (lun_sll=14,lunw=15,mcp=3300000)

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
     & grid     ! Frequency spacing out output spectrum (cm-1)

      parameter (sol=3.0E+08)  ! Speed of Light (m/s)

      write(*,*)' Enter starting and ending frequencies (cm-1):'
      read(*,*) nus, nue
      write(*,*)' Enter frequency point spacing (cm-1):'
      read(*,*) grid
      write(*,*)' Fraction of solar diameter observed '
      write(*,*)' (0=disc-center; 1=disk-integrated)'
      read(*,*) frac
      ncp=1+(nue-nus)/grid  ! number of spectral points
      ifsp=nint(nus/grid)
      if(ncp.gt.mcp) then
         write(*,*) 'ncp,mcp=',ncp,mcp
         stop ' Increase parameter MCP'
      endif

c      osvc=2500.0d0   ! Observer-Sun velocity component (m/s)
      osvc=0.0d0      ! Observer-Sun velocity component (m/s)

      write(*,'(a)')
     & 'Enter the name of solar linelist (e.g. solar_merged.108): '
      read(*,'(a)')solarll

      call get_ggg_environment(gggpath,dl)
      lr = lnbc(gggpath)

      call solar_pts(lun_sll,gggpath(:lr)//dl//'linelist'//dl//solarll,
     & ifsp,grid*(1.d0+osvc/sol),frac,spts,ncp)

      if(nus.le.999.5 .and. nue.le.9999.5) then
         fmt='(a1,i3.3,a1,i4.4,a1,i3.3,a4)'
      elseif(nus.le.9999.5 .and. nue.le.9999.5) then
         fmt='(a1,i4.4,a1,i4.4,a1,i3.3,a4)'
      elseif(nus.le.999.5 .and. nue.le.99999.5) then
         fmt='(a1,i3.3,a1,i5.5,a1,i3.3,a4)'
      elseif(nus.le.9999.5 .and. nue.le.99999.5) then
         fmt='(a1,i4.4,a1,i5.5,a1,i3.3,a4)'
      else
         fmt='(a1,i5.5,a1,i5.5,a1,i3.3,a4)'
      endif
   
      outfile=solarll
      ldot=index(outfile,'.')
      write(outfile(ldot:),fmt) '_',nint(nus),'_',nint(nue),'_',
     & nint(100*frac),'.out'
      open(lunw,file=outfile, status='unknown')
      do i=1,ncp
         write(lunw,'(f8.2,f8.5)') (ifsp+i-1)*grid,spts(i)
      end do
      close(lunw)
      stop
      end
