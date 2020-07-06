c  adjust_solar_linelist.f
c  Program to revise a solar linelist based on fits to measured spectra.
c  The solar lines are assumed to have an absorption shape:
c    h(x,d,w) = Exp2[-x^2/SQRT(d^4+x^2.w^2(1+|x|/(w+d))]
c  The solar pseudo transmittance is then  T = rspf + (1-rspf)*Exp1(-s.h )
c  where s is the line strength and rspf is BB(STmin)/BB(T(v)), the Ratio
c  of the Solar Planck Functions at the T minimum (STmin=4000K) and at
c  unit optical depth ~5800K). This relates to the Minnaert correction.
c  Note that there are two different Exp operators, one converting
c  the absorption to transmittance, and the other describing the
c  (Gaussian) lineshape. To avoid confusion, these are labeled Exp1
c  and Exp2.
c
c  The partial differentials are:
c    dT/ds = -h.(1-rspf).Exp1(-s.h)
c    dT/dh = -s.(1-rspf).Exp1(-s.h)
c    dh/dx =  h.x.(2d^4+x^2.w^2(1+0.5|x|/(w+d)))/den^1.5
c    dh/dd =  h.x^2.(2d^3-0.5|x|.(w.x/(w+d))^2)/den^1.5
c    dh/dw =  h.w*x^4(1+|x|.(w+2d)/(w+d)^2)/den^1.5
c  where den=(d^4+x^2.w^2)
c
c  s.dT/ds =             = -s.(1-rspf).Exp1(-s.h).h
c    dT/dx = dT/dh.dh/dx = -s.(1-rspf).Exp1(-s.h).dh/dx
c    dT/dd = dT/dh.dh/dd = -s.(1-rspf).Exp1(-s.h).dh/dd
c    dT/dw = dT/dh.dh/dw = -s.(1-rspf).Exp1(-s.h).dh/dd
c
c  Let b = -s.(1-rspf).Exp1(-s.h).h
c  Let z= |x|/(w+d)
c  s.dT/ds = b
c    dT/dx = b.x.(2d^4+x^2.w^2(1+0.5z))/den^1.5
c    dT/dd = b.x^2.(2d^3-0.5z^3.(w+d).w^2)/den^1.5
c    dT/dw = b.w.x^4(1+z.(w+2d)/(w+d))/den^1.5
c
c  NFP=7 variables are retrieved for each line (frequency,
c  strength_dc, strength_di, wwidth_dc, wwidth_di, dwidth_dc, dwidth_di). 
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
c Program solves the simultaneous equations:
c      dT/dX . Delta_X = Res
c        s_dc / erats =  rats . s_di / erats
c        w_dc / eratw =  ratw . w_di / eratw
c        d_dc / eratd =  ratw . d_di / eratd

c       (s_dc+del_sdc)/erats =  rats . (s_di+del_sdi)/erats
c
c  Run by typing
c  ~/ggg/bin/adjust_solar_linelist window_runlog.col
 
      implicit none

      include "params.f"

      integer*8 file_size_in_bytes,fsib
      integer*4  lun_col,luns,lunt,lunll,nfp,nfpi,nmp,krank,iline,
     & j,k,ncol,
     & reclen,i1,i2,posnall,nline,mline,nrec,kk,ispec,nspec,mspec,
     & lr,nlhead,lnbc
      parameter (nfp=7,mline=680,    ! 500 is the limit under ifort
     & luns=15,lunll=16,lunt=17,mspec=21)
      integer*4 ip(nfp*mline),mw(mline),np(mspec)
      real*8 nu1(mspec),nu2(mspec),f,tm,tc,graw(mspec),freq(mline),
     & xzo,cont,nus,nue,totsdc,totsdi
      real*4 x,x2,hh,den,bb,gg,pd(mmp+2*nfp*mline,nfp*mline),
     & sdc(mline),sdi(mline),wdc(mline),wdi(mline),
     & ddc(mline),ddi(mline),
     & frqcen,width,
     & sct,rspf,
     & rms,rmsbar,
     & stren,smax,
     & frac,dwid,wwid,ff,
     & frstep,d_max,d_min,w_min,
     & res(mmp+2*nfp*mline),wk(nfp*mline),rnorm,effres,asza,peff,
     & ws,ww,wd,tt,stmin,zero,
     & d2,d4, rats,erats, ratw,eratw, ratd,eratd, zobs,zmin
      parameter (smax=99.99999)
      parameter (stmin=4000.0) ! Solar temperature minimum (~500km above surface)
      parameter (rats=1.20,erats=0.8) ! Line strength DC/DI ratio = 1.15+/-0.8
      parameter (ratw=1.00,eratw=0.8) ! Lorentz width DC/DI ratio = 1.00+/-0.8
      parameter (ratd=0.85,eratd=0.4) ! Doppler width DC/DI ratio = 0.85+/-0.4
      parameter (w_min=0.0)
      parameter (zero=0.0)
      parameter (d_max=9.999)

      character llformat*23,version*40,vstamp*12,
     & sss(mline)*37, inputfile*80
c
      llformat='(i3,f13.6,6f9.5,1x,a37)'

      version= ' adjust_solar_linelist  2019-12-30   GCT'
      vstamp='ASL'//version(25:28)//version(30:31)//version(33:35)
      vstamp(12:12)=char(10)
      write(*,*) version
      write(*,*) vstamp

      rmsbar=0.50  ! assumed % average RMS spectral fit
c
c  Read .col file.
      if (iargc() == 0) then
         lun_col = 5
      elseif (iargc() == 1) then
         call getarg(1, inputfile)
         lun_col = 10
         open(lun_col, file=inputfile, status='old')
      else
         stop 'Usage: $gggpath/bin/adjust_solar_linelist colfile'
      endif
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
      
c  Compute Ratio of Solar Planck Functions (for Minnaert correction)
      if(index(solarll,'minnaert').eq.0) then
         rspf=0.0
      else
         sct=2200+1000*log10(frqcen)  ! Solar Continuum Temperature
         rspf=(exp(1.4388*frqcen/sct)-1)/(exp(1.4388*frqcen/stmin)-1)
      endif

c  Open linelist and read relevent portion into memory.
      lr=lnbc(solarll)
      read(solarll(lr-2:lr),*)reclen
      fsib=file_size_in_bytes(lunll,solarll)
      nrec=int(fsib/reclen,kind(nrec))
      if ( nrec*reclen .ne. fsib ) then
         write(*,*)'Linelist size not divisible by record length',reclen
         write(*,*)solarll,fsib
         stop
      endif

      write(*,*) 'Number of records:',nrec,solarll
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
         if(sdc(iline)*sdi(iline).le.0.0) then
            close(lunll)
            write(*,*)'Freq SDC SDI=',freq(iline),sdc(iline),sdi(iline)
            stop 'Error: Line strengths have opposite sign or are zero'
         endif
         if(abs(sdc(iline)).le.0.0 .and. abs(sdi(iline)).le.0.0) then
            close(lunll)
            stop 'Error: Both line strength are zero!'
         endif
         if(abs(sdc(iline)).le.0.0) sdc(iline)=sdi(iline)/4
         if(abs(sdi(iline)).le.0.0) sdi(iline)=sdc(iline)/4
      end do
c
      kk=0
      nspec=0
      do ispec=1,mspec  !  loop over spectra
22       read(lun_col,*,end=99)specname
         if(specname.eq.'ace-solar-spectrum.txt')   goto 22
         if(specname.eq.'fsunallp.500000.txt.wn')   goto 22
         if(specname.eq.'830619R0.005')             goto 22 ! high airmass
         if(specname.eq.'901218R0.003')             goto 22 ! high airmass
         if(specname.eq.'iz20070510NI.00')          goto 22 ! high airmass
         if(specname.eq.'pa20041111saaaaa.002_003') goto 22 ! high airmass
         if(specname.eq.'pa20041111saaaaa.226_227') goto 22 ! high airmass
         if(specname.eq.'pa20041222saaaaa.019')     goto 22 ! high airmass
         if(specname.eq.'pa20041222saaaaa.020')     goto 22 ! high airmass
         if(specname.eq.'db20070417seccaa.206_207') goto 22 ! high airmass
         if(specname.eq.'db20070629seccaa.008_009') goto 22 ! high airmass
         if(specname.eq.'oc20090121sfddaa.346_349') goto 22 ! high airmass
         sptfile=sptpath(:lnbc(sptpath))//specname
         open(luns,file=sptfile,status='old')
         read(luns,*)nlhead,ncol
         read(luns,*)nu1(ispec),nu2(ispec),np(ispec),
     &   effres,asza,zobs,zmin,rms,peff,frac,xzo
c         write(*,*)nu1(ispec),nu2(ispec),np(ispec),
c     &   effres,asza,zobs,zmin,rms,peff,frac,xzo
         write(*,*) ispec,np(ispec),sptfile
         if(frac.gt.1.) frac=1.
         ff=frac**2
         graw(ispec)=(nu2(ispec)-nu1(ispec))/(np(ispec)-1)
         if(abs(rms).gt. 10.0) rms=0.1*rms**2  ! de-weight really bad fits
         read(luns,*)
         if (kk+np(ispec).gt.mmp) then
            write(*,*)'nmp,mmp=',kk+np(ispec),mmp
            stop 'nmp>mmp'
         endif
         do k=1,np(ispec)  !  loop over points
            kk=kk+1
            read(luns,*)f,tm,tc,cont
            tm=(tm/cont-xzo)/(1-xzo)
            tc=(tc/cont-xzo)/(1-xzo)
            res(kk)=sngl((tm-tc)*tc/(abs(rms)+rmsbar))
c            write(*,*)'f,tm,tc,res(kk)=',f,tm,tc,res(kk)
            do iline=1,nline
               stren=sdc(iline)*(1-ff)+sdi(iline)*ff
               dwid=ddc(iline)*(1-ff)+ddi(iline)*ff
               wwid=wdc(iline)*(1-ff)+wdi(iline)*ff
               d2=dwid**2
               d4=d2**2
               x=sngl(f-freq(iline))
               x2=x**2
               den=d4+x2*wwid**2
               hh=exp(-(x2/sqrt(den)))
c
c  Let b = -s.(1-rspf).Exp1(-s.h).h
c  s.dT/ds =             = b
c    dT/dx = dT/dh.dh/dx = b.x.(d^4+den)/den^1.5
c    dT/dd = dT/dh.dh/dd = b.x^2.2d^3/den^1.5
c    dT/dw = dT/dh.dh/dw = b.w.x^4/den^1.5
               bb=stren*(rspf-1.)*exp(-stren*hh)*hh/(abs(rms)+rmsbar)
               gg=bb/den**1.5
               pd(kk,nfp*(iline-1)+1)=gg*x*(d4+den)          ! d/dx
               pd(kk,nfp*(iline-1)+2)=bb*(1-ff)              ! d/ds1
               pd(kk,nfp*(iline-1)+3)=bb*ff                  ! d/ds2
               pd(kk,nfp*(iline-1)+4)=gg*wwid*x2**2*(1-ff)   ! d/dw1
               pd(kk,nfp*(iline-1)+5)=gg*wwid*x2**2*ff       ! d/dw2
               pd(kk,nfp*(iline-1)+6)=gg*2*dwid**3*x2*(1-ff) ! d/dd1
               pd(kk,nfp*(iline-1)+7)=gg*2*dwid**3*x2*ff     ! d/dd2
            end do   !  iline=1,nline
         end do    ! k=1,np(ispec)
         close(luns)
      end do  ! ispec=1,mspec
      read(lun_col,*,end=99) specname
      write(*,*)'Warning: ispec > mspec = ',ispec,mspec
99    if(iargc() .eq. 1) close(lun_col)
      nspec=ispec-1
      nmp=kk       ! total number of measured spectral points (all used spectra)
c
c  Write out the PD's to file (in XYPLOT format).
c  for first solar line in window only.
      write(*,*) ' nline, nfp = ', nline,nfp
      kk=0
      open(lunt,file='adjust_ll_pd.out',status='unknown')
      write(lunt,*)2,10
      write(lunt,*)' i  f  res  x  sdc  sdi  wdc  wdi  ddc  ddi '
      do ispec=1,nspec
         do k=1,np(ispec)
            kk=kk+1
            write(lunt,'(i2,f11.4,e12.4,7e12.4)')  ispec,
     &      nu1(ispec)+(k-1)*graw(ispec),res(kk),(pd(kk,j),j=1,nfp)
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

c  Add constraint that couples Strengths (DC & DI)
c  Strong lines are de-weighted
         tt=exp(-sngl(sdc(iline)))
         ws=tt/erats
         pd(nmp+nfpi+2,nfpi+2)=ws                      ! set diagonal element
         pd(nmp+nfpi+3,nfpi+2)=-ws                     ! off-diagonal element
c         res(nmp+nfpi+2)=-ws*sngl(sdc(iline)-rats*sdi(iline))
         res(nmp+nfpi+2)=-ws*sngl(1-rats*sdi(iline)/sdc(iline))
         pd(nmp+nfpi+2,nfpi+3)=-rats*ws                ! off-diagonal element
         pd(nmp+nfpi+3,nfpi+3)=rats*ws                 ! set diagonal element
c         res(nmp+nfpi+3)= ws*sngl(sdc(iline)-rats*sdi(iline))
         res(nmp+nfpi+3) = ws*sngl(sdc(iline)/sdi(iline)-rats)

c  Add constraint that couples W-Widths (DC & DI)
         ww=tt/eratw
         pd(nmp+nfpi+4,nfpi+4)=ww                      ! set diagonal element
         pd(nmp+nfpi+5,nfpi+4)=-ww                     ! off-diagonal element
         res(nmp+nfpi+4)=-ww*sngl(wdc(iline)-ratw*wdi(iline))
         pd(nmp+nfpi+4,nfpi+5)=-ratw*ww                ! off-diagonal element
         pd(nmp+nfpi+5,nfpi+5)=ratw*ww                 ! set diagonal element
         res(nmp+nfpi+5)= ww*sngl(wdc(iline)-ratw*wdi(iline))

c  Add constraint that couples D-Widths (DC & DI)
         wd=tt/eratd
         pd(nmp+nfpi+6,nfpi+6)=wd                      ! set diagonal element
         pd(nmp+nfpi+7,nfpi+6)=-wd                     ! off-diagonal element
         res(nmp+nfpi+6)=-wd*sngl(ddc(iline)-ratd*ddi(iline))
         pd(nmp+nfpi+6,nfpi+7)=-ratd*wd                ! off-diagonal element
         pd(nmp+nfpi+7,nfpi+7)=ratd*wd                 ! set diagonal element
         res(nmp+nfpi+7)= wd*sngl(ddc(iline)-ratd*ddi(iline))

      end do  ! iline=1,nline

c  Augment PD with diagonal constraint matrix. Augment RES with zeros.
c  This reduces the step size and tilts it toward the steepest descent direction.
      do iline=1,nline
         nfpi=nfp*(iline-1)
         do j=1,nfp
            call vmov(zero,0,pd(nline*nfp+nmp+1,nfpi+j),1,nfp*nline) ! set column nfpi+j to zero beyond row NMP
            pd(nline*nfp+nmp+nfpi+j,nfpi+j)=var                      ! set diagonal elements to VAR
            res(nline*nfp+nmp+nfpi+j)=0.0
         end do
      end do

c  Solve matrix equation PD.dx=resids
c      write(*,*)' Calling hfti'
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
      frstep=0.5
      totsdc=0.0d0
      totsdi=0.0d0
      do iline=1,nline
         nfpi=nfp*(iline-1)
         d_min=4.3E-07*sngl(freq(iline))*sqrt(stmin/mw(iline))
         write(*,'(f12.4,8f12.6)')freq(iline),(res(nfpi+j),j=1,nfp)
         write(25,'(f12.4,8f12.6)')freq(iline),(res(nfpi+j),j=1,nfp)
         freq(iline)=freq(iline)+dble(frstep*res(nfpi+1))
c         sdc(iline)=sdc(iline)+frstep*dble(res(nfpi+2))
c         sdi(iline)=sdi(iline)+frstep*dble(res(nfpi+3))
         sdc(iline)=sdc(iline)*exp(frstep*res(nfpi+2))
         sdi(iline)=sdi(iline)*exp(frstep*res(nfpi+3))
         if(sdc(iline).gt.smax) then
            write(*,*)'Setting sdc to smax',freq(iline),sdc(iline),smax
            sdc(iline)=smax
         endif
         if(sdi(iline).gt.smax) then
            write(*,*)'Setting sdi to smax',freq(iline),sdi(iline),smax
            sdi(iline)=smax
         endif
         wdc(iline)=max(w_min,wdc(iline)+frstep*res(nfpi+4))
         wdi(iline)=max(w_min,wdi(iline)+frstep*res(nfpi+5))
         ddc(iline)=
     &   max(d_min,min(d_max,ddc(iline)+frstep*res(nfpi+6)))
         ddi(iline)=
     &   max(d_min,min(d_max,ddi(iline)+frstep*res(nfpi+7)))
         write(lunll,llformat,rec=iline+i1) mw(iline),freq(iline),
     &   sdc(iline),sdi(iline), wdc(iline),wdi(iline),
     &   ddc(iline),ddi(iline),sss(iline)(:25)//vstamp
         totsdc=totsdc+sdc(iline)
         totsdi=totsdi+sdi(iline)
      end do
      close(lunll)
      write(*,*)' Sum: sdc,sdi = ',totsdc,totsdi
      stop
      end
