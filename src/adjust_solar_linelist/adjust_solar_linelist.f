c  adjust_solar_linelist.f
c  Program to revise a solar linelist based on fits to measured spectra.
c  The solar lines are assumed to have an absorption shape:
c    k = s.exp[-x^2/SQRT(d^4+x^2.y^2)]
c  The solar pseudo transmittance is then  T = exp[-k]
c  The partial differentials are;
c   dk/ds = k/s = exp[-x^2/SQRT(d^4+x^2.y^2)]
c   dk/dx = x.(2d^4+x^2.y^2)/den
c   dk/dd = -2*x^2.d^3/den
c   dk/dy = -y*x^4/den
c  where den=(d^4+x^2.y^2)^1.5
c
c   dT/ds - dT/dk . dk/ds = -T.exp[-x^2/SQRT(d^4+x^2.y^2)]
c
c  NFP=7 variables are retrieved for each line (frequency,
c  strength_dc, strength_di, lwidth_dc, lwidth_di, dwidth_dc, dwidth_di). 
c
c  There are NSPEC * NPTS spectral measurements.
c
c  There are two types of constraint:
c  1) Relating the disk-center and the disk-integrated parameter values
c  2) A more general  Levenberg-Marquardt type constraint that nudges
c   the step update towards steepest descent rather than Gauss-Newton
c
c  The revised values are written back into the same linelist,
c  overwriting the old values.
c
c  Run by typing
c  ~/ggg/bin/adjust_solar_linelist < window_runlog.col
c
      implicit none
      include "../ggg_const_params.f"
      include "params.f"

      integer*4  lun_col,luns,lunt,lunll,nfp,nfpi,nmp,krank,iline,
     & i,j,k,ncol,
     & reclen,i1,i2,posnall,nline,mline,nrec,kk,ispec,nspec,mspec,
     & lr,fsib,nlhead,lnbc,file_size_in_bytes
      parameter (nfp=7,mline=500,lun_col=5,luns=15,lunll=16,
     & lunt=17,mspec=19)
      integer*4 ip(nfp*mline),mw(mline),np(mspec)
      real*8 nu1(mspec),nu2(mspec),f,tm,tc,rms,frqcen,width,graw(mspec)
      real*8 freq(mline),sdc(mline),sdi(mline),wdc(mline),wdi(mline),
     & ddc(mline),ddi(mline),rc,stren,wwid,dwid,
     & frstep,frac,nus,nue,d_min,w_min
      real*4 x,x2,xx,den,ff,ss,gg,pd(mmp+2*nfp*mline,nfp*mline),
     & res(mmp+2*nfp*mline),wk(nfp*mline),rnorm,
     & scons,sconw,scond,ws,ww,wd,tt,
     & d2,d4, rats, ratw, ratd,zobs,zmin,wx2
      parameter (scons=3.0,sconw=2.0,scond=2.0,
     & rats=1.1,  ! Ratio of the line strengths (DC/DI))
     & ratw=0.96, ! Ratio of the line Lorentz widths (DC/DI))
     & ratd=0.99, ! Ratio of the line Doppler widths (DC/DI)
     & w_min=0.001)
      character llformat*23,
     & sss(mline)*37
c
      write(*,*)' adjust_solar_linelist    2011-02-12   GCT'
      llformat='(i3,f13.6,6f9.5,1x,a37)'
c
c  Read .col file.
      read(lun_col,*)nlhead,ncol
      do j=2,nlhead-7
         read(lun_col,*) 
      end do
      read(lun_col,'(34x,a)')solarll
      read(lun_col,*)
      read(lun_col,'(34x,a)')sptpath
      read(lun_col,*)
      read(lun_col,*)
      read(lun_col,*)frqcen,width
      read(lun_col,*)
      nus=frqcen-width/2
      nue=frqcen+width/2
      
c  Open linelist and read relevent portion into memory.
      lr=lnbc(solarll)
      read(solarll(lr-2:lr),*)reclen
      fsib=file_size_in_bytes(lunll,solarll)
      nrec=fsib/reclen
      if ( nrec*reclen .ne. fsib ) then
         write(*,*)'Linelist size not divisible by record length',reclen
         write(*,*)solarll,fsib
         stop
      endif

      write(*,*) lr,reclen,nrec,solarll
      open(lunll,file=solarll, access='direct',form='formatted',
     & status='old',recl=reclen)
c  Don't fit lines that are centered outside the window
c  They will be better fitted in the previous/next panel
      i1=posnall(lunll,nus,nrec)
      i2=posnall(lunll,nue,nrec)
      nline=i2-i1
      if(nline.gt.mline) then
         write(*,*)'nline,mline=',nline,mline
         stop 'nline.gt.mline'
      endif
      do iline=1,nline
         read(lunll,llformat,rec=iline+i1) mw(iline),freq(iline),
     &   sdc(iline),sdi(iline),
     &   wdc(iline),wdi(iline),
     &   ddc(iline),ddi(iline),sss(iline)
      end do
c
      kk=0
      do ispec=1,mspec  !  loop over spectra
         read(lun_col,*,end=99)specname
         sptfile=sptpath(:lnbc(sptpath))//specname
         write(*,*) ispec,sptfile
         open(luns,file=sptfile,status='old')
         read(luns,*)nlhead,ncol
         read(luns,*)nu1(ispec),nu2(ispec),np(ispec),
     &   zobs,zmin,zmin,rms,frac
         if(frac.gt.1.) frac=1.
         ff=frac**2
         graw(ispec)=(nu2(ispec)-nu1(ispec))/(np(ispec)-1)
         if(abs(rms).gt.9.9) rms=2.0E+34
         read(luns,*)
         if (kk+np(ispec).gt.mmp) then
            write(*,*)'nmp,mmp=',kk+np(ispec),mmp
            stop 'nmp>mmp'
         endif
         do k=1,np(ispec)  !  loop over points
            kk=kk+1
            read(luns,*)f,tm,tc
            res(kk)=(tm-tc)/(abs(rms)+0.02)
            do iline=1,nline
               stren=sdc(iline)*(1-ff)+sdi(iline)*ff
               dwid=ddc(iline)*(1-ff)+ddi(iline)*ff
               wwid=wdc(iline)*(1-ff)+wdi(iline)*ff
               d2=dwid**2
               d4=d2**2
               x=f-freq(iline)
               x2=x**2
               wx2=x2*wwid**2
               den=d4+wx2
               xx=exp(-(x2/sqrt(den)))
c               rc=0.70138-3.8252d-5*freq(k)  ! Minnaert correction.
               rc=0.0
               ss=(rc-1.)*tc*xx/(abs(rms)+0.02)
               gg=ss*stren/den**1.5
               pd(kk,nfp*(iline-1)+1)=x*(d4+den)*gg              ! d/dx
               pd(kk,nfp*(iline-1)+2)=ss*(1-ff)                  ! d/ds1
               pd(kk,nfp*(iline-1)+3)=ss*ff                      ! d/ds2
               pd(kk,nfp*(iline-1)+4)=wwid*x2**2*gg*(1-ff)   ! d/dw1
               pd(kk,nfp*(iline-1)+5)=wwid*x2**2*gg*ff       ! d/dw2
               pd(kk,nfp*(iline-1)+6)=2*dwid**3*x2*gg*(1-ff) ! d/dd1
               pd(kk,nfp*(iline-1)+7)=2*dwid**3*x2*gg*ff     ! d/dd2
            end do   !  iline=1,nline
         end do    ! k=1,np(ispec)
         close(luns)
      end do  ! ispec=1,mspec
      read(lun_col,*,end=99) specname
      write(*,*)'Warning: nspec > mspec'
99    nspec=ispec-1
      nmp=kk       ! total number of measured spectral points (from all spectra)
      close(lun_col)
c
c  Write out the PD's to file (in XYPLOT format).
        write(*,*) ' nline, nfp = ', nline,nfp
        kk=0
        open(lunt,file='adjust_ll_pd.out',status='unknown')
        write(lunt,*)2,10
        write(lunt,*)' i  f  res  x  sdc  sdi  wdc  wdi  ddc  ddi '
        do ispec=1,nspec
        do k=1,np(ispec)
           kk=kk+1
           write(lunt,'(i2,f11.4,f9.4,7e12.4)')  ispec,
     &     nu1(ispec)+(k-1)*graw(ispec),res(kk),(pd(kk,j),j=1,nfp)
        end do
        end do
      close(lunt)

c  Augment PD with off-diagonal constraints that relate the DC and DI values. 
c     First, set the lower part of PD with zeros
      do iline=1,nline
         nfpi=nfp*(iline-1)
         do j=1,nfp
            call vmov(zero,0,pd(nmp+1,nfpi+j),1,nfp*nline) ! set column nfpi+j to zero beyond row NMP
         end do

         tt=exp(-sdc(iline))
c  Add constraint that couples Strengths (DC & DI)
         ws=tt*scons
         pd(nmp+nfpi+2,nfpi+2)=ws                      ! set diagonal element
         pd(nmp+nfpi+3,nfpi+2)=-ws                     ! off-diagonal element
         res(nmp+nfpi+2)=-ws*(sdc(iline)-rats*sdi(iline))
         pd(nmp+nfpi+2,nfpi+3)=-rats*ws                ! off-diagonal element
         pd(nmp+nfpi+3,nfpi+3)=rats*ws                 ! set diagonal element
         res(nmp+nfpi+3)= ws*(sdc(iline)-rats*sdi(iline))

c  Add constraint that couples W-Widths (DC & DI)
         ww=tt*sconw
         pd(nmp+nfpi+4,nfpi+4)=ww                      ! set diagonal element
         pd(nmp+nfpi+5,nfpi+4)=-ww                     ! off-diagonal element
         res(nmp+nfpi+4)=-ww*(wdc(iline)-ratw*wdi(iline))
         pd(nmp+nfpi+4,nfpi+5)=-ratw*ww                ! off-diagonal element
         pd(nmp+nfpi+5,nfpi+5)=ratw*ww                 ! set diagonal element
         res(nmp+nfpi+5)= ww*(wdc(iline)-ratw*wdi(iline))

c  Add constraint that couples D-Widths (DC & DI)
         wd=tt*scond
         pd(nmp+nfpi+6,nfpi+6)=wd                      ! set diagonal element
         pd(nmp+nfpi+7,nfpi+6)=-wd                     ! off-diagonal element
         res(nmp+nfpi+6)=-wd*(ddc(iline)-ratd*ddi(iline))
         pd(nmp+nfpi+6,nfpi+7)=-ratd*wd                ! off-diagonal element
         pd(nmp+nfpi+7,nfpi+7)=ratd*wd                 ! set diagonal element
         res(nmp+nfpi+7)= wd*(ddc(iline)-ratd*ddi(iline))

      end do  ! iline=1,nline

c  Augment PD with diagonal constraint matrix. Augment RES with zeros.
c  This reduces the step size and tilts it toward the steepesst descent direction.
      do iline=1,nline
         nfpi=nfp*(iline-1)
         do j=1,nfp
            call vmov(zero,0,pd(nline*nfp+nmp+1,nfpi+j),1,nfp*nline) ! set column nfpi+j to zero beyond row NMP
            pd(nline*nfp+nmp+nfpi+j,nfpi+j)=var                      ! set diagonal elements to VAR
            res(nline*nfp+nmp+nfpi+j)=0.0
         end do
      end do

c  Solve matrix equation PD.dx=resids
      call shfti(pd,mmp+2*nfp*mline,nmp+2*nfp*nline,nfp*nline,
     & res,nmp+2*nfp*mline,1,tau,krank,rnorm,wk,ip)
      if(krank.lt.nfp*nline) write(6,*)'Rank Deficient:',krank,' /',
     & nfp*nline
c
      write(*,*)'     f          x          sdc        sdi       
     & wdc       wdi         ddc      ddi'
      write(25,*)'     f          x          sdc        sdi       
     & wdc       wdi         ddc      ddi'

c  Update values and write into linelist.
      frstep=0.5d0
      do iline=1,nline
         nfpi=nfp*(iline-1)
         d_min=4.3E-07*freq(iline)*sqrt(4000.0/mw(iline))
         write(*,'(f12.4,8f12.6)')freq(iline),(res(nfpi+j),j=1,nfp)
         write(25,'(f12.4,8f12.6)')freq(iline),(res(nfpi+j),j=1,nfp)
         freq(iline)=freq(iline)+frstep*res(nfpi+1)
         sdc(iline)=sdc(iline)+frstep*dble(res(nfpi+2))
         sdi(iline)=sdi(iline)+frstep*dble(res(nfpi+3))
         wdc(iline)=max(w_min,wdc(iline)+frstep*dble(res(nfpi+4)))
         wdi(iline)=max(w_min,wdi(iline)+frstep*dble(res(nfpi+5)))
         ddc(iline)=
     &        max(d_min,min(d_max,ddc(iline)+frstep*dble(res(nfpi+6))))
         ddi(iline)=
     &        max(d_min,min(d_max,ddi(iline)+frstep*dble(res(nfpi+7))))
         sss(iline)(37:37)=char(10)
         write(lunll,llformat,rec=iline+i1) mw(iline),freq(iline),
     &   sdc(iline),sdi(iline),
     &   wdc(iline),wdi(iline),ddc(iline),ddi(iline),sss(iline)
      end do
      close(lunll)
      stop
      end
