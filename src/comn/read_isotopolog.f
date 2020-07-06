      subroutine read_isotopolog(lunr,iso_fmt,kgas,kiso,gas,speci_id,
     & icode,fia,delta,lnfrd,molewt,ewvb,atc,tdrpf,vibfrq,dgen,nmode,
     & mmode,istat)
c
c  Reads one record of the (already opened) isotopologs.dat file
c
c  Inputs:
c         lunr  I*4   Logical unit number 
c      iso_fmt  C*120 Format of data on isotopologs.dat
c        mmode  I*4   Declared dimension of arrays vibfrq and degen
c
c  Outputs:  
c         kgas  I*4   Gas ID number
c         kiso  I*4   Isotopolog
c          gas  C*8   Gas name (e.g.CH3CN)
c     speci_id  C*24  Full gas name
c          fia  R*4   Fractional Isotopic Abundance
c        delta  R*4   Isotopic fractionation (per mil) in troposphere
c        lnfrd  R*4   Ln(fr) dependence of Isotopic fractionation 
c       molewt  I*4   Molar Mass (g)
c          atc  R*4   Additional Temperature Correction (non-zero for PLL only)
c        tdrpf  R*4   Temperature Dependence of Rotational Partition Function
c vibfrq(nmode) R*4   Fundamental Vibrational Frequencies
c  dgen(nmode)  I*4   Degeneracies
c        nmode  I*4   Number of fundamental Vibrational modes
c

      implicit none
      integer*4  j,
     &   lunr,          ! Logical Unit number
     &   kgas,          ! Gas ID number
     &   kiso,          ! ISotopolog
     &   icode,         ! Isotopolog code
     &   nmode,         ! Number of different vibrational modes (ie 3N-6) of each gas
     &   mmode,         ! Maximum possible value of NMODE
     &   dgen(mmode),   ! degeneracy of vibrational modes
     &   molewt,        ! Molar Mass
     &   istat          ! Status Flag
c
      real*4
     &   ewvb,
     &   atc,           ! Additional Temperature Correction (for PLL)
     &   tdrpf,         ! T-Dependence of Rotational Partition Function
     &   fia,           ! Fractional isotopic abundance
     &   delta,         ! Isotopic fractionation in troposphere
     &   lnfrd,         ! Gradient of Isotopic fractionation wrt ln(f)
     &   vibfrq(mmode)  ! Array of vibrational frequencies
c
      character
     &   iso_fmt*(*),
     &   col1*1,
     &   gas*(*),       ! Gas name
     &   speci_id*(*)   ! Full name
c
c
      col1=':'
      istat=0
c  Skip commented out lines
c      do while (col1.eq.':' .or. col1.eq.';')  
1     read(lunr,iso_fmt,iostat=istat, end=99)
     &col1,kgas,kiso,gas,speci_id,icode,fia,delta,lnfrd,
     &molewt,ewvb,atc,tdrpf,nmode,(vibfrq(j),dgen(j),j=1,nmode)
      if(col1.eq.':' .or. col1.eq.';') go to 1

c      write(*,*)'read_isotop: istat,kgas,kiso,nmode =',istat,kgas,
c     &kiso,nmode
      do j=1,nmode
c         write(*,*)'Vibfrq: kgas,kiso,j,vibfrq,dgen,istat=',
c     &   kgas,kiso,j,vibfrq(j),dgen(j),istat
         if(vibfrq(j).le.0.0 ) then
            write(*,*)'Warning: read_isotopologs: Invalid Vibfrq:'//
     &      'kgas,kiso,j,vibfrq,dgen=',kgas,kiso,j,vibfrq(j),dgen(j)
         endif
      end do
c      write(*,*)kgas,kiso,gas

      if(nmode.gt.mmode) then
         write(*,*)'NMODE= MMODE= ',nmode, mmode
         write(*,*)'(Use isotopologs.dat instead of isotopomers.dat)'
         stop 'READ_ISOTOP: Increase parameter MMODE'
      endif
99    continue
      return
      end
