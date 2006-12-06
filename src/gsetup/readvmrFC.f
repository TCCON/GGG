      subroutine readvmrFC(lunr,vmrpath,z,nlev,pabel,vmr,mgas,modname)
c  Reads and interpolates (possibly with a vertical shift) an initial vmr set
c  (SETNAME) onto the vertical grid defined by Z(NLAY) and places the result
c  into array VMR.

c     Extended from readvmr.f to also read FASCODE format vmr files
c     DG Sept00

      implicit none

      INTEGER*4 lunr,klev,nlev,jcol,ncol,mgas,ngas,mg,nn,nss,k,pos,
     &          minlvl, ninlvl, igas, i, ii, lnbc
c      parameter (mg=134,minlvl=151)    !minlvl <= mg  -DG: not necessary if inlvl and vold not equivalenced.
      parameter (mg=230,minlvl=250)    !DG Jan03
      CHARACTER vmrpath*(*),pabel*(*),dum*16,string*1000,modname*48
      real*4 z(nlev),zold,znew,vold(mg),vnew(mg),vmr(mgas,nlev),fr,
     &vshift,inlvl(minlvl)
c      equivalence (inlvl,vold)
      if(mgas.gt.mg) then
        write(*,*) 'mgas=',mgas
        write(*,*) 'mg=',mg
        stop 'readvmrFC: increase parameter MG'
      endif
      vshift=0.0
c==================================================================
c  Read in names of gases in vmr set
c      write(6,*)'readvmr: mgas,nlev ',mgas,nlev
      open(lunr,file=vmrpath,status='old')
      read(lunr,'(a)')string
c
      if(index(vmrpath,'.vmr').gt.0.or.index(vmrpath,'.set').gt.0) then          !DG Jan03
c         GFIT format
          backspace (lunr)
          read(lunr,*)nn,ncol
          ngas=ncol-1   ! first column is Z
          do k=2,nn
            read(lunr,'(a)')pabel
            pos = index(pabel,'VSHIFT:') 
            if( pos .gt. 0 ) read(pabel(pos+7:),*) vshift
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
            stop 'READVMR: mismatch b/n number of data values & labels'
          endif
          read(string,*)zold,(vold(jcol),jcol=1,ngas)
c          write(*,*)'readvmrFC:',zold
c
          read(lunr,'(a)')string
          call substr(string,dum,1,nss)
          if(nss.ne.ncol) then
            write(*,*) nss,ncol
            stop 'READVMR: mismatch b/n number of data values & labels'
          endif
          read(string,*)znew,(vnew(jcol),jcol=1,ngas)
c          write(*,*)'readvmrFC:',znew
c
          if(z(1).lt.zold+vshift) then
            write(6,*)'Warning: vmrs may not extend low enough'
            write(*,*) 'z(1), zold, vshift', z(1), zold, vshift
          endif
          do klev=1,nlev
c          write(*,*)klev,nlev,z(klev),znew
 11         if(z(klev).gt.znew+vshift) then
              zold=znew
              do jcol=1,ngas
                vold(jcol)=vnew(jcol)
              end do
              read(lunr,'(a)',end=110)string
              call substr(string,dum,1,nss)
              if(nss.ne.ncol) then
                write(*,*) nss,ncol,klev
                stop 'READVMR: mismatch b/n no. of data values & labels'
              endif
              read(string,*)znew,(vnew(jcol),jcol=1,ngas)
              go to 11
            else
c             write(*,*)fr,znew,zold
              fr=(z(klev)-zold-vshift)/(znew-zold)
              do jcol=1,ngas
                vmr(jcol,klev)=fr*vnew(jcol)+(1.-fr)*vold(jcol)
              end do
             endif
          end do
          close(lunr)
c 33       write(6,*)'Warning: vmrs do not extend high enough'
c          stop
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
110   return
      end
