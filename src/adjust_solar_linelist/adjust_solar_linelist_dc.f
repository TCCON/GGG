c  adjust_solar_linelist.f
c  Program to revise a solar linelist based on fits to measured spectra.
c  The solar lines are assumed to have an absorption shape:
c    k = s.exp[-x^2/SQRT(d^4+x^2.y^2)]
c  The solar pseudo transmittance is then  T = exp[-k]
c  The partial differentials are;
c   dk/ds = k/s
c   dk/dx = x.(2d^4+x^2.y^2)/den
c   dk/dd = -2*x^2.d^3/den
c   dk/dy = -y*x^4/den
c  where den=(d^4+x^2.y^2)^1.5
c
c  Four variables are retrieved for each line (strength, frequency,
c  doppler width, Lorentz width). 
c
c  The revised values are written back into the same linelist,
c  overwriting the old values.
c
c  Run by typing
c  ~/ggg/bin/adjust_solar_linelist window_runlog.col
c
      implicit none
      include "../ggg_const_params.f"
      include "params.f"

      integer*4 lunr,luns,lunw,nfp,nmp,krank,i,j,k,ncol,
     & reclen,k1,k2,posnall,nline,mline,nrec,kk,ispec,nspec,mspec,
     & lr,fsib,nlhead,lnbc,file_size_in_bytes
      parameter (nfp=4,mline=600,luns=15,lunw=16,
     & mspec=10)
      integer*4 ip(nfp*mline),mw(mline),np(mspec)
      real*8 nu1(mspec),nu2(mspec),f,tm,tc,rms,frqcen,width,graw(mspec)
      real*8 freq(mline),stren(mline),w_wid(mline),d_wid(mline),rc,
     & frstep,sbhw(mline),eprime(mline),tdpbhw(mline),nus,nue,
     & d_min,w_min
      real*4 x,x2,den,ff,ss,gg,pd(mmp+nfp*mline,nfp*mline),
     & res(mmp+nfp*mline),wk(nfp*mline),rnorm,
     & d2,d4,zobs,zmin,wx2
      parameter (d_min=0.005,w_min=0.0001)
      character llformatr*46,llformatw*54,llformatx*54,
     & sss(mline)*37, inputfile*50
c
      write(*,*)' solar_linelist    2009-11-13   GCT'
      llformatr='(i3,f12.6,e10.3,e10.3,2f5.4,f10.4,f8.4,1x,a37)'
      llformatw='(i3,f12.6,1pe10.3,e10.3,0pf5.4,f5.4,f10.4,f8.4,1x,a37)'
      llformatx='(i3,f12.6,1pe10.3,e10.3,0pf5.3,f5.4,f10.4,f8.4,1x,a37)'
c
c  Read .col file.
      if (iargc() == 0) then
         lunr = 5
      elseif (iargc() == 1) then
         lunr = 10
         call getarg(1, inputfile)
         open(lunr, file=inputfile, status='old')
      else
         stop 'Usage: $gggpath/bin/adjust_solar_linelist colfile'
      endif
      read(lunr,*)nlhead,ncol
      do j=2,nlhead-6
         read(lunr,*) 
      end do
      read(lunr,'(34x,a)')solarll
      read(lunr,*)
      read(lunr,'(34x,a)')sptpath
      read(lunr,*)
      read(lunr,*)frqcen,width
      read(lunr,*)
      nus=frqcen-width/2
      nue=frqcen+width/2
      
c  Open linelist and read relevent portion into memory.
      lr=lnbc(solarll)
      read(solarll(lr-2:lr),*)reclen
      fsib=file_size_in_bytes(lunw,solarll)
      nrec=fsib/reclen
      if ( nrec*reclen .ne. fsib ) then
         write(*,*)'Linelist size not divisible by record length',reclen
         write(*,*)solarll,fsib
         stop
      endif

      write(*,*) lr,reclen,nrec,solarll
      open(lunw,file=solarll, access='direct',form='formatted',
     & status='old',recl=reclen)
c  Don't fit lines that are centered outside the window
c  The will be better fitted in the previous/next panel
      k1=posnall(lunw,nus,nrec)
      k2=posnall(lunw,nue,nrec)
      nline=k2-k1
      if(nline.gt.mline) then
         write(*,*)'nline,mline=',nline,mline
         stop 'nline.gt.mline'
      endif
      do k=1,nline
         read(lunw,llformatr,rec=k+k1) mw(k),freq(k),stren(k),
     &   w_wid(k),d_wid(k),sbhw(k),eprime(k),tdpbhw(k),sss(k)
      end do
c
      i=0
      do ispec=1,mspec  !  loop over spectra
         read(lunr,*,end=99)specname
         sptfile=sptpath(:lnbc(sptpath))//specname
         open(luns,file=sptfile,status='old')
         read(luns,*)nlhead,ncol
         read(luns,*)nu1(ispec),nu2(ispec),np(ispec),zobs,zmin,zmin,rms
         graw(ispec)=(nu2(ispec)-nu1(ispec))/(np(ispec)-1)
         if(abs(rms).gt.9.9) rms=2.0E+34
         read(luns,*)
         if (i+np(ispec).gt.mmp) then
            write(*,*)'nmp,mmp=',i+np(ispec),mmp
            stop 'nmp>mmp'
         endif
         do j=1,np(ispec)  !  loop over points
            i=i+1
            read(luns,*)f,tm,tc
            res(i)=(tm-tc)/rms
            kk=0
            do k=1,nline
               d2=d_wid(k)**2
               d4=d2**2
               x=f-freq(k)
               x2=x**2
               wx2=x2*w_wid(k)**2
               den=d4+wx2
               ff=exp(-(x2/sqrt(den)))
               rc=0.70138-3.8252d-5*freq(k)
               rc=0.0
               ss=(rc-1.)*tc*ff/rms
               gg=ss*stren(k)/den**1.5
               pd(i,1+kk)=x*(d4+den)*gg                      ! d/dx
               pd(i,2+kk)=ss                                 ! d/ds
               pd(i,3+kk)=w_wid(k)*x2**2*gg  ! d/dw
               pd(i,4+kk)=2*d_wid(k)**3*x2*gg       ! d/dd
               kk=kk+nfp
            end do   !  k=1,nline
         end do    ! j=1,np(ispec)
         close(luns)
      end do  ! ispec=1,mspec
      read(lunr,*,end=99)specname
      write(*,*)'Warning: nspec > mspec'
99    nspec=ispec-1
      nmp=i            ! total number of measured spectral points (from all spectra)
      if (iargc() == 1) close(lunr)
c
c  Write out the PD's to file (in XYPLOT format).
         write(*,*) nline,nfp
        i=0
        write(25,*)2,7
        write(25,*)' i f  res  x  s  w  d '
        do ispec=1,nspec
        do j=1,np(ispec)
           i=i+1
           write(25,'(i2,f11.4,f9.4,4e12.4)')  ispec,
     &     nu1(ispec)+(j-1)*graw(ispec),res(i),(pd(i,k),k=1,nfp)
        end do
        end do

      do i=1,nfp*nline
         call vmov(zero,0,pd(1+nmp,i),1,nfp*nline)
         pd(i+nmp,i)=var
         res(i+nmp)=0.0
      end do 
c
c  Solve matrix equation PD.dx=resids
      call shfti(pd,mmp+nfp*mline,nmp+nfp*nline,nfp*nline,
     & res,nmp+nfp*mline,1,tau,krank,rnorm,wk,ip)
      if(krank.lt.nfp*nline) write(6,*)'Rank Deficient:',krank,' /',
     & nfp*nline
c
      write(*,*)'     f          x          s          w         d'
c  Updates values and write into linelist.
      frstep=0.5d0
      kk=0
      do k=1,nline
         write(*,'(f12.4,8f12.6)')freq(k),(res(kk+j),j=1,nfp)
         freq(k)=freq(k)+frstep*res(kk+1)
         stren(k)=stren(k)+frstep*dble(res(kk+2))
         w_wid(k)=max(w_min,w_wid(k)+frstep*dble(res(kk+3)))
         d_wid(k)=max(d_min,min(d_max,d_wid(k)+frstep*dble(res(kk+4))))
         sss(k)(37:37)=char(10)
         if(d_wid(k).le.0.9999) then
            write(lunw,llformatw,rec=k+k1) mw(k),freq(k),stren(k),
     &      w_wid(k),d_wid(k),sbhw(k),eprime(k),tdpbhw(k),sss(k)
         else
            write(lunw,llformatx,rec=k+k1) mw(k),freq(k),stren(k),
     &      w_wid(k),d_wid(k),sbhw(k),eprime(k),tdpbhw(k),sss(k)
         endif
         kk=kk+nfp
      end do
      close(lunw)
      stop
      end
