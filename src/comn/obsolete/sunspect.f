      subroutine sunspect(solarll,nlines,reclen,fzero,grid,frac,sot,ncp)
c
c  Calculates the Solar Absorption spectrum.

c  INPUTS:
c     SOLARLL  C*(*) Path to solar linelist
c     NLINES   I*4   Number of lines in solar linelist (e.g. 17362)
c     RECLEN   I*4   Record Length of SOLARLL (e.g. 101)
c     FZERO    R*8   Frequency of zero'th spectral point (cm-1)
c     GRID     R*8   Spectral Point spacing (cm-1)
c     FRAC     R*8   Fraction of the solar diameter viewed
c     NCP      I*4   Number of points to be calculated
c
c OUTPUTS:
c     SOT(NCP) R*4  SOLAR Optical Thickness spectrum
c
c
c DESCRIPTION:
c  Calculates the solar optical thickness spectrum (SOT) at frequencies
c       Vi = fzero + i * grid       i = 1,NCP
c  and adds it to the contents of array SOT(NCP).
c
c  Taking the exponential of SOT produces the solar spectrum
c  as it would be observed at infinite spectral resolution.
c
c  All solar lines are assumed to have a shape of the form
c
c          SOT = s.exp(-x^2/sqrt(d^4+x^2.y^2))
c  where
c          s is the line-center optical thickness (dimensionless)
c          x is the frequency from line center (cm-1)
c          y is the 1/e folding width (cm-1)
c          d is the Doppler width (cm-1)
c
c  In the doppler limit, i.e. d^2 >> x.y  
c         SOT = s.exp(-(x/d)^2)
c
c  In the far line wing limit, i.e. x.y >> d^2,  
c         SOT = s.exp(-|x/y|)
c
c  So near the line center, the lineshape is Doppler, but in
c  the line wings it decays exponentially (if y>0).
c
c  This choice of lineshape has no physical basis. It just seems
c  to give a reasonable representation is nearly all cases.
c  The only cases in which this lineshape does not give an
c  adequate representation of the absorption are the extremely
c  broad lines of light atmos such as H (atomic hydrogen) or Mg.
c  However, by representing the H absorptions as superpositions
c  of two lines, one narrow and the other broad, adequate results
c  were obtained.
c
c  Molecular absorptions (e.g. CO, OH, NH, CN) tend to have narrow,
c  Doppler lineshapes because they are confined to a relatively
c  narrow layer in the cooler, upper, part of the solar atmosphere.
c  In the hotter depths they are dissociated.
c
c  Atomic transitions, on the other hand, are formed over a much
c  wider range of solar altitudes, and hence temperatures. This
c  gives rise to line shapes whose wings decay in an approximately
c  exponential manner with the distance from line center. The line
c  shape of equation (1) does a reasonable job in both cases.
c
c  This subroutine also makes allowances for the effect of the
c  finite FOV of the observing instrument, which gives rise to:
c  (1) broadening of the solar lines due to the linear variation
c  of the Doppler shift from solar rotation across the solar disk.
c  (2) deepening of the solar lines due to limb darkening.
c  It assumes that an instrument which observes the entire solar
c  disk will observe lines which are, on average, twice the strength
c  of an instrument just observing the center of the disk.
c
c  Note that array SOT is NOT initialized in this subroutine,
c  so the solar optical thickness spectrum is added to whatever
c  is already there. This allows it to be added to the atmospheric
c  optical thickness spectrum and the sum exponentiated together.
c
      implicit none
      integer*4 mw,iline,ncp,kline1,kline2,kv1,kv2,iv,
     & nlines,posnall,reclen
      real*4 sot(ncp),smax,yy
      real*8 fzero,grid,flinwid,sld,srot,frac,xx,x2,d4,y2,margin,
     & ss,rr,pi,acc,freq,w_wid,stren,d_wid,sbhw,eprime,tdpbhw
      character llformat*53,sss*36,solarll*(*)
      parameter (pi=3.14159265d0,acc=0.0001d0,margin=20.,smax=26.0)
      llformat='(i3,f12.6,2e10.3,f5.4,f5.4,f10.4,f8.4,1x,a36)'
c
      if ( ncp .lt. 1 ) stop ' SOLARSPEC: NCP < 1   '
c
      sld=2.d0/(1.d0+sqrt(1.d0-frac**2))   ! Solar limb darkening ??
      sld=1.0

c  Open linelist
      open(19,file=solarll,access='direct',
     & form='formatted',status='old',recl=reclen)
      kline1=posnall(19,fzero-margin,nlines)
      kline2=posnall(19,fzero+grid*ncp+margin,nlines)
      do iline=kline1+1,kline2
         read(19,llformat,rec=iline) mw,freq,stren,w_wid,d_wid,sbhw,
     &   eprime,tdpbhw,sss
         srot=5.E-06*freq*frac    ! broadening due to solar rotation
         d4=(d_wid**2+srot**2)**2  ! Total Gaussian width
         flinwid=sqrt(2*stren*(d_wid+w_wid)/acc)
         kv1=max0(1,int((freq-fzero-flinwid)/grid))
         kv2=min0(ncp,int((freq-fzero+flinwid)/grid))
         y2=(w_wid)**2
         xx=grid*kv1-(freq-fzero)
         ss=sld*stren
         do iv=kv1,kv2
            xx=xx+grid
            x2=xx**2
            rr=x2/sqrt(d4+y2*x2*(1+8.5*abs(xx)))
            rr=x2/sqrt(d4+y2*x2*(1+abs(xx/(w_wid+0.07))))
            yy=ss*exp(-rr)
            sot(iv)=sot(iv)-yy
         end do
      end do
99    close(19)
      return
      end
