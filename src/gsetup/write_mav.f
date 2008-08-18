      subroutine write_mav(z,t,p,d,vmr,nlev,lun_mav,
     & t_cell,p_cell,gas_in_cell,vmrlabel,isofile,mgas)
c
c  writes the model and vmr profiles to LUNV
      implicit none
      integer*4 lun_mav,nspeci,jspeci,nlev,i,j,kspeci,zdiff,
     & i10,mmode,istat,jvmr,kvmr,nvmr,lunt,mgas,gas_in_cell
      parameter (kspeci=230,kvmr=70,mmode=30,lunt=20)
      integer*4 gasindex(kspeci),sagindex(kvmr),
     & kgas,kiso,dumint,ifail,
     & dumintarr(mmode),nmode,speci_cell,nspeci_cell
      real*4 z(nlev+1),t(nlev+1),p(nlev+1),d(nlev+1),vmr(mgas,nlev+1),
     &     delta(kspeci),epsilon(kspeci),fia,dumreal,
     &     dumrealarr(mmode),boltzmann,t_cell,p_cell,zero
      character isofile*(*),
     &     shortname(kspeci)*(8), fullname(kspeci)*(11),
     &     vmrlabel*(*),vmrgases(kvmr)*10,speci_id*24

      boltzmann=1.38066E-23
      zero=0.0
      data sagindex/kvmr*0/

      open(19,file='write_mav.rpt',status='unknown')
      speci_cell=0
      nspeci_cell=0
c     Index gas names in vmrlabel to gas names in "isotopologs.dat"
c      write(6,*)' Opening ',isofile
      open(lunt,file=isofile,status='old')
      do jspeci=1,kspeci
         call read_isotop(lunt,kgas,kiso,shortname(jspeci),speci_id,
     &   fia,delta(jspeci),epsilon(jspeci),dumint,dumreal,
     &   dumrealarr,dumintarr,nmode,mmode,istat)
         if(istat.ne.0) goto 77
         if(kgas.eq.gas_in_cell) then
            speci_cell=jspeci
            nspeci_cell=nspeci_cell+1
         endif
         call lowercase(shortname(jspeci))
         fullname(jspeci)=char(kiso+ichar('0'))//shortname(jspeci)
         gasindex(jspeci)=0
c         write(*,*) jspeci,kgas,shortname(jspeci), speci_id,
c     &        fia,delta(jspeci),epsilon(jspeci),
c     & dumint, dumreal, nmode
      end do
      write(6,*) 'write_mav: Increase parameter KSPECI'
      stop
 77   close(lunt)
      nspeci=jspeci-1
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
         write(*,*)' See write_mav.rpt for details'
         stop 'GSETUP Error'
      endif
c     
c     Output model information (SUNRUN.MAV)
      write(lun_mav,*) 2, nspeci+4, nlev+1
      write(lun_mav,'(a38,230a11)')
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
c First line of .mav file contains cell information
          write(lun_mav,'(2f7.2,2(1pe11.3),230e11.3)')
     $    -9.9,t_cell,p_cell/1013.25,0.0001*p_cell/boltzmann/t_cell,
     &    (zero,j=1,speci_cell-nspeci_cell),
     &    (1.0,j=1,nspeci_cell),
     &    (zero,j=speci_cell+1,nspeci)
      do i=1,nlev
         if (z(i).le.10.0) then
c     Write tropospheric fractions
            write(lun_mav,'(2f7.2,2(1pe11.3),144e11.3)')
     $           z(i),t(i),p(i),d(i),(vmr(gasindex(j),i)*
     $           (1.0+delta(j)/1000.0),j=1,nspeci)
         else
c     Write stratospheric enrichments          
            write(lun_mav,'(2f7.2,2(1pe11.3),144e11.3)')
     $           z(i),t(i),p(i),d(i),((1.0+delta(j)/1000.0)*
     &           vmr(gasindex(j),i)**(1.+epsilon(j)/1000.0)/
     $           vmr(gasindex(j),i10)**(epsilon(j)/1000.0),
     &           j=1,nspeci)
         end if
      end do      ! i=1,nlev
      return
      end
