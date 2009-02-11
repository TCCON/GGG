      subroutine solar_spectrum(solarll,fzero,grid,frac,spts,ncp)
c
c  Calculates a pseudo-transmittance spectrum of the sun.
c  The spectrum can be disk center or disk-integrated
c  dependening on the selected linelist (SOLARLL).

c  INPUTS:
c     SOLARLL  C*(*) Path to solar linelist
c     FZERO    R*8   Frequency of zero'th spectral point (cm-1)
c     GRID     R*8   Spectral Point spacing (cm-1)
c     FRAC     R*8   Fraction of the solar diameter viewed
c     NCP      I*4   Number of points to be calculated
c
c OUTPUTS:
c     STS(NCP) R*4   Solar Transmittance Spectrum
c
c
c DESCRIPTION:
c  Calculates the solar Pseudo-Transmittance Spectrum (STS) at frequencies
c       Vi = fzero + i * grid       i = 1,NCP
c  where fzero and grid are specified by the user in the calling program.
c
c  It is recommended that the spectral point spacing not exceed
c  the doppler widths of the narrowest solar features of interest. 
c  For solar CO, these widths are about 0.02 cm-1 at 2100 cm-1
c  and 0.04 cm-1 at 4200 cm-1
c  This spectrum is computed as it would be observed at infinite spectral
c  resolution and must therefore eventually be convolved with an ILS.
c

c BASIS FOR LINESHAPE FORMULATION
c  Molecular absorptions (e.g. CO, OH, NH, CN) tend to have narrow,
c  Doppler lineshapes because they are confined to a relatively
c  narrow layer in the cooler, upper, part of the solar atmosphere.
c  In the hotter depths molecules are dissociated.
c
c  Atoms, on the other hand, are stable to much higher temperatures.
c  Atomic absorptions can occur at greater depths in the sun where
c  the Temperature and pressure are much larger, giving wider
c  absorption features. The net result of absorption at different
c  depths inside the sun is a line shapes whose wings decay in an
c  approximately exponential manner with the distance from line center.
c
c  The line shape of below does a reasonable job for both atoms
c and molecules.
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
c  Note that array STS is NOT initialized in this subroutine,
c  so the solar optical thickness spectrum is added to whatever
c  is already there. This allows it to be added to the atmospheric
c  optical thickness spectrum and the sum exponentiated together.
c
      implicit none
      integer*4 mw,iline,ncp,kline1,kline2,kv1,kv2,iv,i,lr,lnbc,
     & nlines,posnall,reclen,fsib,file_size_in_bytes
      real*4 spts(ncp),zero
      real*8 fzero,grid,flinwid,sld,srot,frac,xx,x2,d4,y2,margin,
     & ss,rr,ff,pi,acc,freq,w_wid,stren,d_wid,sbhw,eprime,tdpbhw,rc
      character llformat*53,sss*36,solarll*(*)
      parameter (pi=3.14159265d0,acc=0.0001d0,margin=20.)
      llformat='(i3,f12.6,2e10.3,f5.4,f5.4,f10.4,f8.4,1x,a36)'
c
c      write(*,*)'SOLAR_SPEC',grid,fzero,ncp
      if ( ncp .lt. 1 ) stop ' SOLARSPEC: NCP < 1   '
c
      sld=1.0-0.05*frac**2  ! solar limb darkening
      sld=1.0               ! solar limb darkening
      zero=0.0
      if(frac.gt.1.0) then
          write(*,*) 'Warning: solar_spectrum: frac > 1',frac
          frac=1.
      endif

c  Determine length of records in solar linelist (RECLEN).
c  This is encoded into the last 3 characters of its name, e.g. 101
      lr=lnbc(solarll)
      read(solarll(lr-2:lr),*)reclen

c Determine the toal size of the solar linelist (FSIB)
c and divide this by RECLEN to determine the number of lines (NLINES).
c Check that NLINES is an integer.
      fsib=file_size_in_bytes(19,solarll)
      nlines=fsib/reclen
      if ( nlines*reclen .ne. fsib ) then
         write(*,*)'Linelist size not divisible by record length',reclen
         write(*,*)solarll,fsib
         stop
      endif

c  Initialize array SPTS to zero
      do iv=1,ncp
          spts(iv)=0.0
      end do

c  Open solar linelist and read lines between fzero and fzero+grid*nCP
      open(19,file=solarll,access='direct',
     & form='formatted',status='old',recl=reclen)
      kline1=posnall(19,fzero-margin,nlines)
      kline2=posnall(19,fzero+grid*ncp+margin,nlines)
c      write(*,*) kline1, kline2
      do iline=kline1+1,kline2
         read(19,llformat,rec=iline) mw,freq,stren,w_wid,d_wid,sbhw,
     &   eprime,tdpbhw,sss
c         write(*,*)iline,freq,stren
         srot=5.E-06*freq*sqrt(frac)    ! broadening due to solar rotation
         d4=(d_wid**2+srot**2)**2  ! Total Gaussian width
         flinwid=sqrt(abs(2*stren*(d_wid+w_wid)/acc))
         kv1=max0(1,int((freq-fzero-flinwid)/grid))
         kv2=min0(ncp,int((freq-fzero+flinwid)/grid))
         y2=(w_wid)**2
         ss=sld*stren
         ff=fzero-freq
         do iv=kv1,kv2
            xx=ff+iv*grid
            x2=xx**2
            rr=x2/sqrt(d4+y2*x2*(1+abs(xx/(w_wid+0.07))))
            spts(iv)=spts(iv)-ss*exp(-rr)
         end do
      end do
99    close(19)
c
c  Convert optical thickness into an apparent transmittance
c  and then apply Minnaert correction (not currently enabled).
      freq=fzero
      do i=1,ncp
        freq=freq+grid
        spts(i)=exp(spts(i))
c        rc=0.70138-3.8252d-5*freq
c        spts(i)=rc+(1.-rc)*spts(i)
      end do 
      return
      end
