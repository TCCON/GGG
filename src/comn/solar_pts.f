      subroutine solar_pts(lunr,solarll,ifsp,grid,frac,spts,ncp)
c
c  Calculates a pseudo-transmittance spectrum of the sun.
c  The spectrum can be disk center or disk-integrated
c  dependening on the selected linelist (SOLARLL).

c  INPUTS:
c     LUNR     I*4   Logical Unit Number to be used
c     SOLARLL  C*(*) Path to solar linelist
c     FZERO    R*8   Frequency of zero'th spectral point (cm-1)
c     GRID     R*8   Spectral Point spacing (cm-1)
c     FRAC     R*8   Fraction of the solar diameter viewed
c     NCP      I*4   Number of points to be calculated
c
c OUTPUTS:
c     SPTS(NCP) R*4   Solar Pseudo-Transmittance Spectrum
c
c
c DESCRIPTION:
c  Calculates the solar Pseudo-Transmittance Spectrum (SPTS) at
c  frequencies
c       V(i) = grid * (ifsp+i-1)       i = 1,NCP
c  where ifsp and grid are specified by the user in the calling program.
c
c  It is recommended that the spectral point spacing not exceed
c  the doppler widths of the narrowest solar features of interest. 
c  For solar CO, these widths are about 0.02 cm-1 at 2100 cm-1
c  and 0.04 cm-1 at 4200 cm-1. Typically, a much-narrower spectral
c  point spacing is used to adequately sample telluric absorptions.
c
c  The solar spectrum is computed as it would be observed at infinite
c  spectral resolution and must therefore eventually be convolved
c  with an ILS.
c
c
c BASIS FOR LINESHAPE FORMULATION
c  Molecular absorptions (e.g. CO, OH, NH, CN) tend to have narrow,
c  Doppler lineshapes because they are confined to a relatively
c  narrow layer in the cooler, upper, part of the solar photosphere.
c  In the hotter depths molecules are dissociated.
c
c  Atoms, on the other hand, are stable to much higher temperatures.
c  Atomic absorptions can occur at greater depths in the sun where
c  the temperature and pressure are much larger, giving wider
c  absorption features. The net result of absorption at different
c  depths inside the sun is "cusp-shaped" lines whose wings decay
c  in an approximately exponential manner with the distance from
c  line center.
c
c  The line shape used here does a reasonable job for both atoms
c  and molecules.
c
c          Absorption = s.exp(-x^2/sqrt(d^4+x^2.y^2))
c  where
c          s is the line-center optical thickness (dimensionless)
c          x is the frequency from line center (cm-1)
c          y is the 1/e folding width (cm-1)
c          d is the Doppler width (cm-1)
c
c  In the doppler limit, i.e. d^2 >> x.y  
c         Absorption = s.exp(-(x/d)^2)
c
c  In the far line wing limit, i.e. x.y >> d^2,  
c         Absorption = s.exp(-|x/y|)
c
c  So near the line center, the lineshape is Doppler, but in
c  the line wings it decays exponentially (if y>0).
c
c  This choice of lineshape has no physical basis. It just seems
c  to give a reasonable representation is nearly all cases.
c  The only exceptions are the extremely broad lines of light
c  atoms such as H (atomic hydrogen) or Mg. By representing
c  these cases as superpositions of two lines, however,
c  adequate results are obtained.
c
c OTHER NOTES
c  This subroutine also makes allowances for the effect of the
c  finite FOV of the observing instrument, which gives rise to:
c  broadening of the solar lines due to the linear variation
c  of the Doppler shift from solar rotation across the solar disk.
c
      implicit none
      integer*4 mw,iline,ncp,kline1,kline2,kv1,kv2,iv,i,lr,mflag,
     & lnbc,nlines,posnall,reclen,lunr,ifsp
      integer*8 fsib,file_size_in_bytes 
      real*4 spts(ncp),zero,aa,rspf,stmin,sct,sfreq
      real*8
     & grid,flinwid,srot,frac,dd,xx,x2,d4,y2,margin,
     & sdc,sdi,wdc,wdi,ddc,ddi,
     & rr,ff,acc,freq,w_wid,stren,d_wid
      character llformat*16,solarll*(*)
      parameter (acc=0.00001d0,margin=90.,stmin=4000.0)
      llformat='(i3,f13.6,6f9.5)'
c
      if(index(solarll,'minnaert').eq.0) then
         mflag=0
      else
         mflag=1
      endif

      if ( ncp .lt. 1 ) stop ' SOLARSPEC: NCP < 1   '
c
      zero=0.0
      if(frac.gt.1.0) then
         write(*,*) 'Warning: solar_spectrum: frac > 1',frac
         frac=1.
      endif
      ff=frac**2

c  Determine length of records in solar linelist (RECLEN).
c  This is encoded into the last 3 characters of its name, e.g. 101
      lr=lnbc(solarll)
      read(solarll(lr-2:lr),*)reclen

c  Determine the total size of the solar linelist (FSIB)
c  and divide this by RECLEN to find the number of lines (NLINES).
c  Check that NLINES is an integer.
      fsib=file_size_in_bytes(lunr,solarll)
      nlines=int(fsib/reclen,kind(reclen))
      if ( nlines*reclen .ne. fsib ) then
         write(*,*)'Linelist size not divisible by record length',reclen
         write(*,*)solarll(:lr),fsib
         stop
      endif
      if( reclen.ne.108) stop 'using wrong solar linelist'

c  Initialize array SPTS to zero
      do iv=1,ncp
         spts(iv)=0.0
      end do

c  Open solar linelist and read lines between
      open(lunr,file=solarll,access='direct',
     & form='formatted',status='old',recl=reclen)
      kline1=posnall(lunr,grid*ifsp-margin,nlines)
      kline2=posnall(lunr,grid*(ifsp-1+ncp)+margin,nlines)
c      write(*,*) kline1, kline2
      do iline=kline1+1,kline2
         read(lunr,llformat,rec=iline) mw,freq,sdc,sdi,wdc,wdi,ddc,ddi
         stren=(1-ff)*sdc+ff*sdi
         w_wid=(1-ff)*wdc+ff*wdi
         d_wid=(1-ff)*ddc+ff*ddi
c         write(*,*)iline,freq,stren
c         srot=5.E-06*freq*sqrt(frac)    ! broadening due to solar rotation
         aa=0.538
         srot=3.95E-06*freq*frac/sqrt(aa+(1-aa)*frac**2.5)    ! broadening due to solar rotation
         d4=(d_wid**2+srot**2)**2  ! Total Gaussian width
         flinwid=dsqrt(dabs(2.d0*stren*(d_wid+w_wid)/acc)) !  Effective line width
         kv1=max0(1,int((freq-flinwid)/grid)-ifsp+1)
         kv2=min0(ncp,int((freq+flinwid)/grid)-ifsp+1)
         y2=(w_wid)**2
         dd=w_wid+4*mflag*d_wid+0.07d0*(1-mflag) !  GCT 20130207
         do iv=kv1,kv2
            xx=grid*(ifsp-1+iv)-freq
            x2=xx**2
c            rr=x2/sqrt(d4+y2*x2*(1+abs(xx/(w_wid+0.07))))
            rr=x2/dsqrt(d4+y2*x2*(1.d0+dabs(xx/dd))) ! GCT 20130207
            spts(iv)=spts(iv)+sngl(stren*exp(-rr))
         end do
      end do
      close(lunr)

c  Convert optical thickness into an apparent transmittance
c
      if(index(solarll,'minnaert').gt.0) then
c  Apply Minnaert correction using Ratio of Solar Planck Functions (RSPF).
         freq=grid*ifsp
         sfreq=sngl(freq)
         do i=1,ncp
            sct=2200+1000*log10(sfreq)  ! Solar Continuum Temperature
            rspf=(exp(1.4388*sfreq/sct)-1)/(exp(1.4388*sfreq/stmin)-1)
            spts(i)=rspf+(1.-rspf)*exp(-spts(i))
            freq=freq+grid
            sfreq=sngl(freq)
         end do 
      else
         do i=1,ncp
            spts(i)=exp(-spts(i))
         end do
      endif
      return
      end
