      subroutine read_refvmrs(lunr,vmrpath,nlev,z,mgas,modname,vmrlabel,
     & refvmr,ngas,strend,gradlat,seacycle,
     & reflat_vmr,refdate_vmr,refztrop_vmr)

c  Reads and interpolates an initial vmr set onto the vertical grid
c  defined by Z(NLAY) and places the result into array REFVMR.

c     Extended from readvmr.f to also read FASCODE format vmr files
c     DG Sept00
c
c  Inputs:
c     lunr              I*4    Logical Unit Number
c     vmrpath           C**    Path to vmr file
c     nlev              I*4    Number of levels
c     z(nlev)           R*4    Altitudes of levels (km)
c     mgas              I*4    Declared first dimension of VMR matrix
c     modname           C*(*)  Name of model (used only by DG code)
c
c Outputs:
c     vmrlabel          C**    Column labels of vmr file
c     refvmr(mgas,nlev) R*4    Gas VMRs
c     ngas              I*4    Actual number of gases in .vmr file
c     reflat_vmr        R*8    latitude of reference vmrs
c     refdate_vmr       R*8    Date of reference vmrs
c     refztrop_vmr      R*8    Tropopause Altitudes from .vmr file
c     gradlat(mgas)     R*4    Latitude gradients (fractional Eq-to-Pole)


      implicit none

      INTEGER*4 lunr,klev,nlev,jcol,ncol,mgas,ngas,mg,nlheader,nss,
     & k,pos,minlvl,ninlvl,igas,i,ii
c      parameter (mg=134,minlvl=151)    !minlvl <= mg  -DG: not necessary if inlvl and vold not equivalenced.
      parameter (mg=230,minlvl=250)    !DG Jan03
      CHARACTER vmrpath*(*),vmrlabel*(*),dum*16,string*1000,modname*(*)
      real*4 z(nlev),zold,znew,vold(mg),vnew(mg),
     & strend(mgas),gradlat(mgas),seacycle(mgas),
     & refvmr(mgas,nlev),fr,
     & inlvl(minlvl),zero
      real*8 twopi,reflat_vmr,refdate_vmr,refztrop_vmr
c      equivalence (inlvl,vold)

      twopi=8*datan(1.0d0)

      if(mg.lt.mgas) then
         write(*,*) 'mgas=',mgas
         write(*,*) 'mg=',mg
         stop 'read_vmr_fc: increase parameter MG'
      endif
      zero=0.0
c==================================================================
c  Read in names of gases in vmr set
c      write(6,*)'readvmr: mgas,nlev ',mgas,nlev
c      write(6,*)'lunr, lnblnk, vmrpath= ',lunr,lnblnk(vmrpath),vmrpath
      open(lunr,file=vmrpath,status='old')
      read(lunr,'(a)')string
c
      if(index(vmrpath,'vmr').gt.0.or.index(vmrpath,'.set').gt.0)then !DG Jan03
c        GFIT format
         backspace (lunr)
         read(lunr,*)nlheader,ncol
         ngas=ncol-1   ! first column is Z
         if(ngas.gt.mgas) stop 'ngas >  mgas'
         refztrop_vmr= -999.0 ! km   default tropopause altitude
         do k=2,nlheader
            read(lunr,'(a)')vmrlabel
            pos = index(vmrlabel,'ZTROP_VMR:')
            if( pos .gt. 0 ) read(vmrlabel(pos+10:),*) refztrop_vmr
            pos = index(vmrlabel,'LAT_VMR:')
            if( pos .gt. 0 ) read(vmrlabel(pos+8:),*) reflat_vmr
            pos = index(vmrlabel,'DATE_VMR:')
            if( pos .gt. 0 ) read(vmrlabel(pos+9:),*) refdate_vmr
            pos = index(vmrlabel,'ZTROP:')
            if( pos .gt. 0 ) read(vmrlabel(pos+6:),*) refztrop_vmr
            pos = index(vmrlabel,'LATIT:')
            if( pos .gt. 0 ) read(vmrlabel(pos+6:),*) reflat_vmr
            pos = index(vmrlabel,'YEAR:')
            if( pos .gt. 0 ) read(vmrlabel(pos+5:),*) refdate_vmr
            pos = index(vmrlabel,'GRADLAT:')
            if( pos .gt. 0 ) read(vmrlabel(pos+8:),*) (gradlat(igas),
     &       igas=1,ngas)
            pos = index(vmrlabel,'STREND:')
            if( pos .gt. 0 ) read(vmrlabel(pos+7:),*) (strend(igas),
     &       igas=1,ngas)
            pos = index(vmrlabel,'SEACYCLE:')
            if( pos .gt. 0 ) read(vmrlabel(pos+9:),*) (seacycle(igas),
     &       igas=1,ngas)
         end do
c
c     Check that the number of column labels (NSS) matches NCOL
         call substr(vmrlabel,dum,1,nss)
         if(nss.ne.ncol) then
            write(*,*)vmrlabel
            write(*,*) ncol,nss
            stop 'READ_REFVMR: mismatch: # of columns/labels differ'
         endif
c
c     Read VMRs level by level testing to make sure that there are NCOL values
         read(lunr,'(a)')string
         call substr(string,dum,1,nss)
         if(nss.ne.ncol) then
            write(*,*) nss,ncol
            stop 'read_refvmr: mismatch: number of data values/labels'
         endif
         read(string,*)zold,(vold(jcol),jcol=1,ngas)
c
         if(nlev.eq.1) then  ! It's a lab spectrum
            znew=1.E+38
c            call vmov(zero,0,vnew,1,ngas)
            do igas=1,ngas
               vnew(igas)=0.0
            end do
         else
            read(lunr,'(a)')string
            call substr(string,dum,1,nss)
            if(nss.ne.ncol) then
               write(*,*) nss,ncol
               stop 'read_refvmr: mismatch: # of data values/labels'
            endif
            read(string,*)znew,(vnew(jcol),jcol=1,ngas)
         endif
c
         if(z(1).lt.zold) then
            write(6,*)'Warning: vmrs may not extend low enough'
            write(*,*) 'z(1), zold = ', z(1), zold
         endif

         do klev=1,nlev
c            write(*,*)klev,nlev,z(klev),znew
            do while (znew.lt.z(klev))
               zold=znew
               do jcol=1,ngas
                  vold(jcol)=vnew(jcol)
               end do
               read(lunr,'(a)',end=110)string
               call substr(string,dum,1,nss)
               if(nss.ne.ncol) then
                  write(*,*) string
                  write(*,*) nss,ncol,klev
                  stop 'READVMR: mismatch:  # of vmr values & labels'
               endif
               read(string,*)znew,(vnew(jcol),jcol=1,ngas)
            end do  !  while (znew.lt.z(klev))
            fr=(z(klev)-zold)/(znew-zold)
            do jcol=1,ngas
               refvmr(jcol,klev)=fr*vnew(jcol)+(1.-fr)*vold(jcol)
c             if(jcol.eq.70)write(*,*)klev,z(klev),refvmr(jcol,klev)
            end do
         end do   ! klev=1,nlev
         close(lunr)
c33       write(6,*)'Warning: vmrs do not extend high enough'
c         stop
      elseif(index(vmrpath,'.ref').gt.0)then                       !DG jan03
c         FASCODE input format
c         vmr file must contain same number of levels as zpt model
c         levels are first read from zpt model file
         close(unit=lunr)
         open(unit=lunr,file=modname,status='old')
         read(lunr,*)ninlvl, ninlvl 
         read(lunr,*)
         read(lunr,*)(inlvl(ninlvl+1-i),i=1,ninlvl)
         close(lunr)
         open(unit=lunr,file=vmrpath,status='old')
c         vmrlabel='Height '
         do igas=1,mgas
            read(lunr,*,end=110)jcol,dum
            if(jcol.ne.igas)then
               write(*,*)'Species numbers do not match: ',jcol, igas
               go to 110
            endif
c            string=vmrlabel(1:lnblnk(vmrlabel))
c            vmrlabel=string(1:lnblnk(string))//' '//dum
            read(lunr,*)(vnew(ninlvl+1-i),i=1,ninlvl)
            zold=inlvl(1)
            znew=inlvl(2)
            ii=2
            if(z(1).lt.zold)
     &      write(6,*)'Warning: vmrs may not extend low enough'
            do klev=1,nlev
111            if(z(klev).gt.znew) then
                  zold=znew
                  ii=ii+1
                  znew=inlvl(ii)
                  go to 111
               else
                  fr=(z(klev)-zold)/(znew-zold)
                  refvmr(igas,klev)=fr*vnew(ii)+(1.-fr)*vnew(ii-1)
               endif
            end do
         enddo
      else
         write(*,*)'VMR format not recognised'
         write(*,*)'VMR name: ',vmrpath
         stop
      endif

110   close(lunr)
      return
      end
