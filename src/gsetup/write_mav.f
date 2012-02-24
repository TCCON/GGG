      subroutine write_mav(z,t,p,d,vmr,nlev,lunw_mav,lncell,
     & t_cell,p_cell,vmr_cell,gas_in_cell,vmrlabel,isofile,mgas)
c
c  writes the model and vmr profiles to LUNV
      implicit none
      include "../ggg_const_params.f"
      include "../ggg_int_params.f"

      integer*4 lunw_mav,nspeci,jspeci,nlev,i,j,kspeci,zdiff,
     & i10,
     & istat,jvmr,kvmr,nvmr,lunr_iso,mgas,
     & mcell,jcell,lncell,gas_in_cell(lncell)
      parameter (kspeci=230,kvmr=75,lunr_iso=57,mcell=8)
      integer*4 gasindex(kspeci),sagindex(kvmr),
     & kgas,kiso,dumint,ifail,
     & dumintarr(mvmode),nmode,
     & speci_cell(mcell),nspeci_cell(mcell)
      real*4 z(nlev+1),t(nlev+1),p(nlev+1),d(nlev+1),vmr(mgas,nlev+1),
     &  delta(kspeci),epsilon(kspeci),fia,dumreal,
     &  dumrealarr(mvmode),boltzmann,
     &  t_cell(lncell),p_cell(lncell),vmr_cell(lncell)
      character isofile*(*),
     &     shortname(kspeci)*(8), fullname(kspeci)*(11),
     &     vmrlabel*(*),vmrgases(kvmr)*10,speci_id*24

      boltzmann=1.38066E-23
      data sagindex/kvmr*0/

      if(lncell.gt.mcell) stop 'write_mav: lncell > mcell'
      open(19,file='write_mav.rpt',status='unknown')
      do jcell=1,lncell
      speci_cell(jcell)=0
      nspeci_cell(jcell)=0
c     Index gas names in vmrlabel to gas names in "isotopologs.dat"
c      write(6,*)' Opening ',isofile
      open(lunr_iso,file=isofile,status='old')
      do jspeci=1,kspeci
         call read_isotop(lunr_iso,kgas,kiso,shortname(jspeci),speci_id,
     &   fia,delta(jspeci),epsilon(jspeci),dumint,dumreal,
     &   dumrealarr,dumintarr,nmode,mvmode,istat)
         if(istat.ne.0) goto 77
         if(kgas.eq.gas_in_cell(jcell)) then
            speci_cell(jcell)=jspeci
            nspeci_cell(jcell)=nspeci_cell(jcell)+1
         endif
         call lowercase(shortname(jspeci))
         fullname(jspeci)=char(kiso+ichar('0'))//shortname(jspeci)
         gasindex(jspeci)=0
c         write(*,*) jspeci,kgas,shortname(jspeci), speci_id,
c     &        fia,delta(jspeci),epsilon(jspeci),
c     &        dumint, dumreal, nmode
      end do
      write(6,*) 'write_mav: Increase parameter KSPECI'
      stop
 77   close(lunr_iso)
      nspeci=jspeci-1
      end do       !  jcell=1,lncell
c      write(*,*)'nspeci=',nspeci

      call lowercase(vmrlabel)
      call substr(vmrlabel,vmrgases,kvmr,nvmr)
      nvmr=nvmr-1               !  The first column of the vmr file is Altitude.
      if(nvmr.gt.mgas) stop ' GSETUP: NVMR > MGAS'
c     
c  Identify which species are HCl (for cell)
c  Check whether all the gas names in "isotopologs.dat" were also
c  found in the .vmr file.
      ifail=0
      call clistindex(nspeci,fullname,nvmr,vmrgases(2),gasindex)
      call clistindex(nspeci,shortname,nvmr,vmrgases(2),gasindex)
      do jspeci=1,nspeci
c        write(*,*)ifail,jspeci,gasindex(jspeci),vmr(gasindex(jspeci),12)
         if(gasindex(jspeci).le.0) then
            write(19,*) ' write_mav: Warning: '//
     &      shortname(jspeci)//' not found in .vmr file '
            ifail=1
         endif
      end do
c
c  Check whether all the gas names in the .vmr file were also
c  found in "isotopologs.dat"
      call clistindex(nvmr,vmrgases(2),nspeci,fullname,sagindex)
      call clistindex(nvmr,vmrgases(2),nspeci,shortname,sagindex)
      do jvmr=1,nvmr
          if(sagindex(jvmr).le.0) write(19,*) ' write_mav: Warning: '//
     &    vmrgases(1+jvmr)//' not found in isotopologs.dat file: '
      end do
      close(19)
      if(ifail.eq.1) then
         write(*,*)' Error: Gases missing from .vmr file'
         write(*,*)' See fort.19 for details'
         stop 'GSETUP Error'
      endif
c     
c     Output model information (SUNRUN.MAV)
c      write(lunw_mav,'(3i4)') 2, nspeci+4, nlev+lncell
      write(lunw_mav,'(a,230a11)')
     & ' Height  Temp   Pres       Density    ',
     $     (fullname(jspeci),jspeci=1,nspeci)

c  Define tropospheric vmr
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
c          write(*,*)jcell,speci_cell(jcell),nspeci_cell(jcell)
          write(lunw_mav,'(2f7.2,2(1pe11.3),230e11.3)')
     $    -11+1.1*float(jcell),t_cell(jcell),p_cell(jcell)/1013.25,
     &    0.0001*p_cell(jcell)/boltzmann/t_cell(jcell),
     &    (zero,j=1,speci_cell(jcell)-nspeci_cell(jcell)),
     &    (vmr_cell(jcell),j=1,nspeci_cell(jcell)),
     &    (zero,j=speci_cell(jcell)+1,nspeci)
      end do      ! jcell=1,lncell

c Last NLEV lines of .mav file contain atmospheric info
      do i=1,nlev
         if (z(i).le.10.0) then
c     Write tropospheric fractions
            write(lunw_mav,'(2f7.2,2(1pe11.3),230e11.3)')
     $           z(i),t(i),p(i),d(i),(vmr(gasindex(j),i)*
     $           (1.0+delta(j)/1000.0),j=1,nspeci)
         else
c  Write stratospheric enrichments          
c  VMR'(i) = VMR(i)*(1+delta/1000)*[VMR(i)/VMR(itrop)]**(1+epsilon/1000)
            write(lunw_mav,'(2f7.2,2(1pe11.3),230e11.3)')
     $           z(i),t(i),p(i),d(i),((1.0+delta(j)/1000.0)*
     &           vmr(gasindex(j),i)**(1.+epsilon(j)/1000.0)/
     $           vmr(gasindex(j),i10)**(epsilon(j)/1000.0),
     &           j=1,nspeci)
         endif
      end do      ! i=1,nlev
      return
      end
