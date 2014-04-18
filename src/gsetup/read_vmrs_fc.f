      subroutine read_vmrs_fc(lunr,vmrpath,z,nlev,pabel,ztrop_mod,
     & vmr,mgas,modname,tlat,idoy,vmr_found)
c  Reads and interpolates (possibly with a vertical stretch) an initial vmr set
c  (SETNAME) onto the vertical grid defined by Z(NLAY) and places the result
c  into array VMR.

c     Extended from readvmr.f to also read FASCODE format vmr files
c     DG Sept00

      implicit none
c      include "../ggg_const_params.f"

      INTEGER*4 lunr,klev,nlev,jcol,ncol,mgas,ngas,mg,nlheader,nss,
     & k,pos,idoy, minlvl, ninlvl, igas, i, ii, lnbc
c      parameter (mg=134,minlvl=151)    !minlvl <= mg  -DG: not necessary if inlvl and vold not equivalenced.
      parameter (mg=230,minlvl=250)    !DG Jan03
      CHARACTER vmrpath*(*),pabel*(*),dum*16,string*1000,modname*48
      logical*4 vmr_found
      real*4 z(nlev),zold,znew,vold(mg),vnew(mg),vmr(mgas,nlev),fr,
     & inlvl(minlvl),
     & tlat,stuclat
      real*8 zeff,ztrop_mod,ztrop_vmr,year_vmr,latit_vmr,pert,twopi
c      equivalence (inlvl,vold)

      twopi=2*dpi

      if(mgas.gt.mg) then
         write(*,*) 'mgas=',mgas
         write(*,*) 'mg=',mg
         stop 'readvmrFC: increase parameter MG'
      endif
c==================================================================
c  Read in names of gases in vmr set
c      write(*,*)'read_vmr_fc: mgas,nlev ',mgas,nlev
c      write(*,*)' lunr, vmrpath = ',lunr,vmrpath
      open(lunr,file=vmrpath,status='old')
       read(lunr,*)nlheader,ncol
c      read(lunr,'(a)')string
c
      if(index(vmrpath,'vmr').gt.0.or.index(vmrpath,'.set').gt.0)then !DG Jan03
c         GFIT format
          backspace (lunr)
          read(lunr,*)nlheader,ncol
          ngas=ncol-1   ! first column is Z
          ztrop_vmr=ztrop_mod  ! km   default tropopause altitude
          do k=2,nlheader
            read(lunr,'(a)')pabel
            pos = index(pabel,'ZTROP:')
            if( pos .gt. 0 ) read(pabel(pos+6:),*) ztrop_vmr
            pos = index(pabel,'YEAR:')
            if( pos .gt. 0 ) read(pabel(pos+5:),*) year_vmr
            pos = index(pabel,'LATIT:')
            if( pos .gt. 0 ) read(pabel(pos+6:),*) latit_vmr
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
c          write(*,*)'readvmrFC:',zold
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

c  stuclat = Stratospheric Tropical Uplift Center Latitude.
c  This is partly driven by solar heating and partly by ocean temperatures.
c  Assumed to vary sinusoidally between 8S on Jan 30 and 8N on Jul 30.
c  This is 3x smaller than the latitude excursion of the sun and six weeks delayed
          stuclat=-8*dcos(twopi*(idoy-30)/365.25)  ! deg.
          stuclat=0.0  ! Keep it simple to start with
          if (vmr_found) then
              pert=0.0
          else
              pert=exp(-((tlat-stuclat)/25.0)**4)
          endif
c          write(*,*)' stuclat=',stuclat, idoy,pert
c  ZEFF is the effective altitude of the vmrs to be used for level KLEV.
c  If ZTROP_mod = ZTROP_vmr, then ZEFF = Z(KLEV) the true altitude.
c  If ZTROP_mod < ZTROP_vmr, then ZEFF > Z(KLEV) the true altitude,
          do klev=1,nlev
c            write(*,*)klev,nlev,z(klev),znew
             if(vmr_found) then
                zeff=z(klev)
             elseif(z(klev).lt.ztrop_mod) then
                zeff=z(klev)*ztrop_vmr/ztrop_mod  ! troposphere
             else
                zeff=z(klev)+(ztrop_vmr-ztrop_mod)*
     &          exp(-(z(klev)-ztrop_mod)/10.0)       ! stratosphere
                zeff=zeff-pert*3.5*ztrop_mod*(z(klev)/ztrop_mod-1)**2
     &          *exp(-(z(klev)-ztrop_mod)/9.0)
c                if(klev.ge.150) write(*,*) klev,z(klev),zeff
             endif
             if(zeff.gt.z(nlev)) zeff=z(nlev)
 11          if(zeff.gt.znew) then
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
              go to 11
            else
              fr=(zeff-zold)/(znew-zold)
              do jcol=1,ngas
                vmr(jcol,klev)=fr*vnew(jcol)+(1.-fr)*vold(jcol)
              end do
             endif
          end do
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
                  vmr(igas,klev)=fr*vnew(ii)+(1.-fr)*vnew(ii-1)
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
