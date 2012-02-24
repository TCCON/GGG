      subroutine read_vmrs_only(lunr,vmrpath,nlev,z,mgas,modname,pabel,
     & refvmr,ngas,reflat_vmr,refdate_vmr,refztrop_vmr)

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
c     pabel             C**    Column labels of vmr file
c     refvmr(mgas,nlev) R*4    Gas VMRs
c     ngas              I*4    Actual number of gases in .vmr file
c     reflat_vmr        R*8    latitude of reference vmrs
c     refdate_vmr       R*8    Date of reference vmrs
c     refztrop_vmr      R*8    Tropopause Altitudes from .vmr file


      implicit none
      include "../ggg_const_params.f"

      INTEGER*4 lunr,klev,nlev,jcol,ncol,mgas,ngas,mg,nlheader,nss,
     & k,pos,minlvl, ninlvl, igas, i, ii, lnbc
c      parameter (mg=134,minlvl=151)    !minlvl <= mg  -DG: not necessary if inlvl and vold not equivalenced.
      parameter (mg=230,minlvl=250)    !DG Jan03
      CHARACTER vmrpath*(*),pabel*(*),dum*16,string*1000,modname*48
      real*4 z(nlev),zold,znew,vold(mg),vnew(mg),refvmr(mgas,nlev),fr,
     &inlvl(minlvl)
      real*8 twopi,reflat_vmr,refdate_vmr,refztrop_vmr
c      equivalence (inlvl,vold)

      twopi=2*dpi

      if(mg.lt.mgas) then
         write(*,*) 'mgas=',mgas
         write(*,*) 'mg=',mg
         stop 'read_vmr_fc: increase parameter MG'
      endif
c==================================================================
c  Read in names of gases in vmr set
c      write(6,*)'readvmr: mgas,nlev ',mgas,nlev
      open(lunr,file=vmrpath,status='old')
      read(lunr,'(a)')string
c
      if(index(vmrpath,'vmr').gt.0.or.index(vmrpath,'.set').gt.0)then !DG Jan03
c         GFIT format
          backspace (lunr)
          read(lunr,*)nlheader,ncol
          ngas=ncol-1   ! first column is Z
          refztrop_vmr= -999.0 ! km   default tropopause altitude
          do k=2,nlheader
            read(lunr,'(a)')pabel
            pos = index(pabel,'ZTROP_VMR:')
            if( pos .gt. 0 ) read(pabel(pos+10:),*) refztrop_vmr
            pos = index(pabel,'LAT_VMR:')
            if( pos .gt. 0 ) read(pabel(pos+8:),*) reflat_vmr
            pos = index(pabel,'DATE_VMR:')
            if( pos .gt. 0 ) read(pabel(pos+9:),*) refdate_vmr
          end do
c
c     Check that the number of column labels (NSS) matches NCOL
          call substr(pabel,dum,1,nss)
          if(nss.ne.ncol) then
            write(*,*)pabel
            write(*,*) ncol,nss
            stop 'READVMR: mismatch between number of columns & labels'
          endif
c
c     Read VMRs level by level testing to make sure that there are NCOL values
          read(lunr,'(a)')string
          call substr(string,dum,1,nss)
          if(nss.ne.ncol) then
             write(*,*) nss,ncol
             stop 'READVMR: mismatch: number of data values/labels'
          endif
          read(string,*)zold,(vold(jcol),jcol=1,ngas)
c
          if(nlev.eq.1) then  ! It's a lab spectrum
             znew=1.E+38
             call vmov(zero,0,vnew,1,ngas)
          else
             read(lunr,'(a)')string
             call substr(string,dum,1,nss)
             if(nss.ne.ncol) then
                write(*,*) nss,ncol
                stop 'READVMR: mismatch: number of data values/labels'
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
                  write(*,*) nss,ncol,klev
                  stop 'READVMR: mismatch b/n # of data values & labels'
                endif
             read(string,*)znew,(vnew(jcol),jcol=1,ngas)
             end do  !  while (znew.lt.z(klev))
             fr=(z(klev)-zold)/(znew-zold)
             do jcol=1,ngas
                refvmr(jcol,klev)=fr*vnew(jcol)+(1.-fr)*vold(jcol)
c                if(jcol.eq.4) write(*,*) klev,zold,znew,z(klev),
c     &          vold(jcol),vnew(jcol),refvmr(jcol,klev)
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
          pabel='Height '
          do igas=1,mgas
              read(lunr,*,end=110)jcol,dum
              if(jcol.ne.igas)then
                write(*,*)'Species numbers do not match: ',jcol, igas
                return
              endif
              string=pabel(1:lnbc(pabel))
              pabel=string(1:lnbc(string))//' '//dum
              read(lunr,*)(vnew(ninlvl+1-i),i=1,ninlvl)
              zold=inlvl(1)
              znew=inlvl(2)
              ii=2
              if(z(1).lt.zold)
     &        write(6,*)'Warning: vmrs may not extend low enough'
              do klev=1,nlev
 111              if(z(klev).gt.znew) then
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

110   return
      end
