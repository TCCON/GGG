      subroutine hitran_to_atmos_gas_numbering(kgas,kiso)
c  Converts HITRAN gas/isotope numbers to ATMOS numbers
c
c  Inputs:
c      KGAS   I*4  HITRAN Gas number
c      KISO   I*4  HITRAN Isotope Number
c
c  Outputs:
c      KGAS   I*4  ATMOS Gas number
c      KISO   I*4  ATMOS Isotope Number
c
c  Notes:
c     if HITRAN gas doesn't exist in ATMOS (e.g., Atomic O, NO+, HOBr)
c     then  KGAS = 0 is returned.
c
c     Input arguments are overwritten by outputs

      implicit none
      integer*4 kgas,kiso,mgashit
      parameter (mgashit=49)
      integer*4 llindex(mgashit)
      save llindex
      data llindex/
c   ATMOS #    HITRAN #  Gas
     &  1,   !   1      H2O
     &  2,   !   2      CO2
     &  3,   !   3      O3
     &  4,   !   4      N2O
     &  5,   !   5      CO
     &  6,   !   6      CH4
     &  7,   !   7      O2
     &  8,   !   8      N2O
     &  9,   !   9      SO2
     & 10,   !  10      NO2
     & 11,   !  11      NH3
     & 12,   !  12      HNO3
     & 13,   !  13      OH
     & 14,   !  14      HF
     & 15,   !  15      HCl
     & 16,   !  16      HBr
     & 17,   !  17      HI
     & 18,   !  18      ClO
     & 19,   !  19      OCS
     & 20,   !  20      H2CO
     & 21,   !  21      HOCl
     & 41,   !  22      N2
     & 28,   !  23      HCN
     & 30,   !  24      CH3Cl
     & 23,   !  25      H2O2
     & 40,   !  26      C2H2
     & 38,   !  27      C2H6
     & 55,   !  28      PH3
     & 36,   !  29      COF2
     & 50,   !  30      SF6
     & 47,   !  31      H2S
     & 46,   !  32      HCOOH
     & 22,   !  33      HO2
     &  0,   !  34      O
     & 27,   !  35      ClONO2
     &  0,   !  36      NO+
     &  0,   !  37      HOBr
     & 39,   !  38      C2H4
     & 56,   !  39      CH3OH
     & 44,   !  40      CH3Br
     & 59,   !  41      CH3CN
     & 31,   !  42      CF4
     &  0,   !  43      C4H2
     &  0,   !  44      HC3N
     &  0,   !  45      H2
     &  0,   !  46      CS
     &  0,   !  47      SO3
     &  0,   !  48      ??
     &  0/   !  49      ??
c
c  H2O is a bit messy because ATMOS decided
c  to designate HDO as a separate molecule.
      if(kgas.le.0) then
         stop 'hitran_to_atmos_gas_numbering: kgas.le.0'
      elseif(kgas.gt.mgashit) then
         write(*,*) 'kgas,mgashit=',kgas,mgashit
         stop 'hitran_to_atmos_gas_numbering: kgas.gt.mgashit'
      elseif(kgas.eq.1 .and. kiso.ge.7) then
         kgas=71
         kiso=1
      elseif(kgas.eq.1 .and. kiso.ge.4) then
         kgas=49
         kiso=kiso-3
      else
         kgas=llindex(kgas)
      endif

      return
      end
