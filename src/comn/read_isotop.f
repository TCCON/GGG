      subroutine read_isotop(lun,kgas,kiso,gas,speci_id,fia,
     & delta,epsilon,molewt,tdrpf,vibfrq,dgen,nmode,mmode,istat)
c
c  Reads one record of the (already opened) isotopologs.dat file
c
c  Inputs:
c          lun  I*4   Logical unit number
c        mmode  I*4   Declared dimension of arrays vibfrq and degen
c
c  Outputs:  
c         kgas  I*4   Gas ID number
c         kiso  I*4   ISotopolog
c          gas  C*8   Gas name (e.g.CH3CN)
c     speci_id  C*24  Full gas name
c          fia  R*4   fractional Isotopic Abundance
c       molewt  I*4   Molar Mass (g)
c        tdrpf  R*4   Temperature Dependence of Rotational Partition Function
c vibfrq(nmode) R*4   Fundamental Vibrational Frequencies
c  dgen(nmode)  R*4   Degeneracies
c        nmode  I*4   Number of fundamental Vibrational modes
c        istat  I*4   Status Flag (0=success; 1=failure)
c

      implicit none
      integer*4  j,
     &   lun,      ! Logical Unit number
     &   kgas,     !  Gas ID number
     &   kiso,     !  ISotopolog
     &   nmode,    ! Number of different vibrational modes (ie 3N-6) of each gas
     &   mmode,    ! Maximum possible value of NMODE
     &   dgen(mmode),  ! degeneracy of vibrational modes
     &   molewt,    ! Molar Mass
     &   istat     ! Status Flag
c
      real*4
     &   tdrpf,        ! T-Dependence of Rotational Partition Function
     &   fia,delta,epsilon,  ! Fractional Isotopic Abundance
     &   vibfrq(mmode) ! Array of vibrational frequencies
c
      character
     &   col1*1,
     &   gas*(*),        ! Gas name
     &   speci_id*(*)    ! Full name
c
c
      col1=':'
      istat=0
c  Skip commented out lines
      do while ((col1.eq.':' .or. col1.eq.';') .and. istat.eq.0) 
      read(lun,
     & '(a1,2i2,1x,a8,a24,f10.8,f8.0,f10.0,1x,i3,f5.2,i3,30(f6.0,i2))',
     & iostat=istat) col1,kgas,kiso,gas,speci_id,fia,delta,epsilon,
     & molewt,tdrpf,nmode,(vibfrq(j),dgen(j),j=1,nmode)
      if(nmode.gt.mmode) then
         write(*,*)'NMODE= MMODE= ',nmode, mmode
         write(*,*)'(Use isotopologs.dat instead of isotopomers.dat)'
         stop 'READ_ISOTOP: Increase parameter MMODE'
      endif
      end do
      return
      end
