c  average_results.f
c  Averages results from different windows of the same gas.
c
c  Arrays YOBS and YERR are physically 1-D for maximum flexibility,
c  but conceptually they are effectively 2-D arrays YOBS(nrow,nwin)
c  and are treated as such in the averaging_with XXX_bias subroutines
      implicit none
      integer idot,irow,j,jj,jtg,k,ktg,lnbc,fnbc,fbc,i1,
     & npp,lrl,totnit,ntc,mlabel,mval,jval,navg,iav,
     & lr,lunm,lunr,lunc,lunv,lunw,luno,mwin,nwin,iwin,
     & ngas,kgas,iwas,
     & mrow,nauxcol,nlhead,nmiss,iyrwas,doywas,mgas,
     & ncol,kcol,jcol,icol,cwas,jav,
     & nrow,nss,nit,mit,iyr,doy,istat,nfound,locnaux,loc,lwas
      parameter (lunr=12)      ! input file (.xsw)
      parameter (lunw=15)      ! output file (.xav)
      parameter (luno=16)      ! outlier file (.xsw)
      parameter (mgas=80)      ! Max number of gases
      parameter (mwin=600)     ! Total number of columns/windows
      parameter (mrow=160000)  ! Max number of output records/spectra
      parameter (mval=3200000) ! Max number of values (NROW * NCOL)
      parameter (nauxcol=19)   ! Number of auxiliary parameters/columns
      parameter (mlabel=18000) ! Max Number of column lable characters

      integer avindx(mgas+1)
      character ans*1,apf*2,cdum*20,avglabel*1200,
     & gfit_version*48,gsetup_version*48,col_string*500,
     & csformat*80,collabel*(mlabel),swfile*80,avfile*80,col1*1,
     & spname_rl*21,specname_col*21,runlog*72,sign(mrow)*1,
     & spectrum(mrow)*20,tabel*80,clab(2*mwin+nauxcol)*17,
     & colfile*40,collate_version*48,ar_version*48,
     & window(mwin)*10,spname_rlwas*21,avlabel(mgas+1)*8

      real*8 airmass,asza,cl,tilt,zlo,fcen,width,zobs,rmin,rmax,
     & fqshift,graw,obslat,obslon,opd,ovcol,rmsfit,
     & r8was,r8year,r8ydiff,year(mrow),
     & trms,lasf,wavtkr,sia,sis,aipl

      real*4
     & ymiss,rew(mrow),cew(mwin),tew,error_sigma,
     & yaux(nauxcol,mrow),
     & yobs(mval),yerr(mval),
     & ybar(mrow),eybar(mrow),
     & bias(mwin),ebias(mwin)

      real*8 zpdtim,tout,pout,hout,tins,pins,hins,fovi,fovo,amal,
     & snr,zenoff,zoff,sg,zpdwas,max_delta_t,delta_t,zmin
      integer bytepw,ifirst,ilast,possp

      ar_version=' average_results  version 1.0.0  2008-10-01  GCT'

      write(*,'(a)')
     & 'Enter name of .xsw file whose contents are to be averaged'
      read(*,'(a)')swfile
      lr=lnbc(swfile)
      avfile=swfile(:lr-1)//'av'

c  Read the entire contents of the .xsw disk file
      open(lunr,file=swfile,status='old')
      open(lunw,file=avfile,status='unknown')
      read(lunr,'(i2,i4,i6)') nlhead,ncol,nrow
      nwin=(ncol-nauxcol)/2
      read(lunr,'(a)') collate_version
      read(lunr,'(a)') gfit_version
      read(lunr,'(a)') gsetup_version
      read(lunr,'(8x,1016(e12.4))') (ymiss,j=1,ncol)
      read(lunr,'(a)') collabel
      do irow=1,mrow
         read(lunr,'(a1,f13.8,18f13.5,800(e12.4))',end=99)
     &   sign(irow),year(irow),(yaux(k,irow),k=2,nauxcol),
     $   (yobs(irow+nrow*(k-1)),yerr(irow+nrow*(k-1)),k=1,nwin)
      end do  !  irow=1,mrow
99    close(lunr)
      if(nrow.ne.irow-1) stop 'NROW mismatch'
c
c
      open(lunt,file='wincomp.'//ans//'av',status='unknown')
      write(lunt,'(2i3)') 2,4
      write(lunt,'(a)')' Window   Mean_Col   Std_Err   Chi2/N'

      open(luno,file='outliers',status='unknown')
      call substr(collabel,clab,2*mwin+nauxcol,nss)
      if(nss.ne.ncol) stop 'NSS .NE. NCOL'
      write(*,*)nrow,nwin,nss
      locnaux=index(clab(nauxcol+1),'_')
      locnaux=index(collabel,clab(nauxcol+1)(:locnaux-1))-1
      cwas=nauxcol-1
      icol=nauxcol+1
      lwas=1
      kgas=1
      do iwin=1,nwin
         loc=index(clab(icol),'_')
         if(loc.eq.0) loc=index(clab(icol),'-')
         if(clab(cwas)(:lwas-1).ne.clab(icol)(:loc-1)) then  ! new gas
            avindx(kgas)=iwin
            avlabel(kgas)=clab(icol)(:loc-1)
            kgas=kgas+1  ! index of window/column to write averages
            if(kgas.gt.mgas) stop 'kgas.gt.mgas'
         endif
         cwas=icol
         lwas=loc
         icol=icol+2
      end do   !  do iwin=2,nwin
      ngas=kgas-1
      avindx(kgas)=iwin
      avlabel(kgas)=clab(icol)(:loc-1)

      do kgas=1,ngas
         navg=avindx(kgas+1)-avindx(kgas)
         write(*,*)avindx(kgas),avindx(kgas+1)-1,navg,' '//avlabel(kgas)
         jj=1+nrow*(avindx(kgas)-1)
         call average_with_mul_bias (ymiss,nrow,navg,yobs(jj),yerr(jj),
     &   ybar,eybar,bias,ebias,rew,cew,tew)
c
c  Report outliers
         do irow=1,nrow
c            do jcol=avindx(kgas),avindx(kgas+1)-1
            jcol=avindx(kgas)
            do jav=1,navg
               jj=irow+nrow*(jcol-1)
               if(yerr(jj).ne.ymiss) then
               error_sigma=(yobs(jj)-ybar(irow)*bias(jav))/tew/yerr(jj)
               if(abs(error_sigma).gt.5.) then
              write(luno,'(3i5,a12,f9.5,a22,i6,a11,a10)')kgas,irow,jcol,
     &         ' Deviation =',
     &         error_sigma,' sigma for spectrum # ',
     &         nint(yaux(4,irow)),' in window ',clab(nauxcol+2*jcol-1)
c               write(luno,*)yobs(jj),yerr(jj),ybar(irow),bias(jav)
               endif
               endif
               jcol=jcol+1
            end do
         end do
c
c  Move the averages values back into the array, overwriting the earlier stuff
         call vmov( ybar,1,yobs(1+nrow*(kgas-1)),1,nrow)
         call vmov(eybar,1,yerr(1+nrow*(kgas-1)),1,nrow)
c
         do iav=1,navg
            write(32,'(a10,3f10.5)')
     &      clab(nauxcol+2*(avindx(kgas)+iav-1)-1),
     &      bias(iav),ebias(iav),cew(iav)
         end do
         write(32,*)
      end do   ! do kgas=1,ngas
c
c      loc=index(clab(nauxcol+1),'_')
c      i1=index(collabel,clab(nauxcol+1)(:loc-1))
c      avglabel=collabel(:i1-1)
c      icol=nauxcol+1
c      cwas=icol
c      lwas=index(clab(cwas),'_')
c      if(lwas.eq.0) lwas=index(clab(cwas),'-')
c      kgas=1
c      navg=1
c      iwas=1
c      do iwin=2,nwin
c         icol=icol+2
c         loc=index(clab(icol),'_')
c         if(loc.eq.0) loc=index(clab(icol),'-')
cc         write(*,*)iwin,' '//clab(icol)
c         if(clab(cwas)(:lwas-1).ne.clab(icol)(:loc-1)) then  ! new gas
c            avglabel=avglabel(:i1)//
c     &      clab(cwas)(:lwas-1)//'  '//clab(cwas)(:lwas-1)//'_error'
c            i1=i1+26
c            write(*,*)iwas,iwin-1,navg,' '//clab(cwas)(:lwas-1)
c            call average_with_mul_bias (ymiss,nrow,navg,
c     &      yobs(1+nrow*(iwas-1)),yerr(1+nrow*(iwas-1)),
c     &      ybar,eybar,bias,ebias,rew,cew,tew)
cc
cc  Look for outliers
c            do irow=1,nrow
c              do jcol=iwas,iwas+navg-1
c                if(yerr(irow+nrow*(jcol-1)).ne.ymiss) then
c                error_sigma=(yobs(irow+nrow*(jcol-1))-ybar(irow)*
c     &          bias(jcol))/tew/yerr(irow+nrow*(jcol-1))
c                if(abs(error_sigma).gt.5.)
c     &          write(luno,'(a12,f9.5,a22,i6,a11,a10)')' Deviation =',
c     &          error_sigma,' sigma for spectrum # ',
c     &          nint(yaux(4,irow)),' in window ',clab(nauxcol+2*jcol-1)
c                endif
c              end do
c            end do
c            call vmov(ybar,1,yobs(1+nrow*(kgas-1)),1,nrow)
c            call vmov(eybar,1,yerr(1+nrow*(kgas-1)),1,nrow)
c            do iav=1,navg
c               write(33,'(a10,3f10.5)') clab(nauxcol+2*(iwas+iav-1)-1),
c     &         bias(iav),ebias(iav),cew(iav)
c            end do
c            write(33,*)
c            navg=1
c            kgas=kgas+1  ! index of window/column to write averages
c            iwas=iwin    ! First window of current gas
c         else        ! Another window of the current gas
c            navg=navg+1  ! number of windows in current average
c         endif
c         cwas=icol
c
c         lwas=loc
c      end do   !  do iwin=2,nwin
c      write(*,*)iwas,iwin-1,navg,' '//clab(cwas)(:lwas-1)
c      avglabel=avglabel(:i1)//
c     &clab(cwas)(:lwas-1)//'  '//clab(cwas)(:lwas-1)//'_error'
c      call average_with_mul_bias (ymiss,nrow,navg,
c     & yobs(1+nrow*(iwas-1)),yerr(1+nrow*(iwas-1)),
c     & ybar,eybar,bias,ebias,rew,cew,tew)
c      call vmov(ybar,1,yobs(1+nrow*(kgas-1)),1,nrow)
c      call vmov(eybar,1,yerr(1+nrow*(kgas-1)),1,nrow)
c      do iav=1,navg
c         write(33,'(a10,3f10.5)')clab(nauxcol+2*(iwas+iav-1)-1),
c     &   bias(iav),ebias(iav),cew(iav)
c      end do
c

c  Write averaged values to file
      open (lunw,file=avfile,status='unknown')
      write(lunw,'(i2,i4,i6)') nlhead+1,nauxcol+2*ngas,nrow
      write(lunw,'(a)') ar_version
      write(lunw,'(a)') collate_version
      write(lunw,'(a)') gfit_version
      write(lunw,'(a)') gsetup_version
      write(lunw,'(8x,219(1pe12.4))') (ymiss,j=1,nauxcol+2*ngas)
      write(lunw,'(a175,60(a8,2x,a14))') collabel(:locnaux),
     & (avlabel(kgas)(:lnbc(avlabel(kgas))),
     &  avlabel(kgas)(:lnbc(avlabel(kgas)))//'_error',kgas=1,ngas)
      do irow=1,nrow
         write(lunw,'(a1,f13.8,18f13.5,200(1pe12.4))')
     &   sign(irow),(yaux(k,irow),k=1,nauxcol),
     &   (yobs(irow+nrow*(k-1)),yerr(irow+nrow*(k-1)),k=1,ngas)
      end do
      close(lunw)
      close(luno)
      stop
      end

      include 'average_with_mul_bias.f'
      include 'average_with_add_bias.f'
      include '/home/toon/ggg/src/comn/substr.f'
      include '/home/toon/ggg/src/comn/lnbc.f'
      include '/home/toon/ggg/src/comn/fnbc.f'
      include '/home/toon/ggg/src/comn/fbc.f'
      include '/home/toon/ggg/src/comn/vsubs.f'
