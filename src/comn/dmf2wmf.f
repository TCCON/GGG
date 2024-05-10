      subroutine dmf2wmf(mgas,ngas,nlev,ip_h2o,vmr)
c
c  Converts an arrary of Dry Mole Fractions into Wet Mole Fractions.
c  Assumes that the water Vapor DMFs are located in vmr(ip_h2o,*)
c  Inputs:
c      MGAS   I*4  Declared first dimension of VMR
c      NGAS   I*4  Number of gases actually used.
c      NLEV   I*4  Number of different gases.
c      IP_H2O I*4  Column number of H2O
c      VMR(MGAS,NLEV)  R*4  Dry Mole Fractions
c
c  Outputs:
c      VMR(MGAS,NLEV)  R*4  Wet Mole Fractions
c
c  Theory
c      Uses the equation
c       WMF(kgas,ilev) = DMF(kgas,ilev) / (1+DMF(ip_h2o,ilev))

      integer*4 mgas,ngas,kgas,nlev,ilev,ip_h2o
      real*4  vmr(mgas,nlev),dmf_h2o

      if( ip_h2o .lt. 1) stop ' IP_H2O < 1'
      if( ip_h2o .gt. ngas) stop ' IP_H2O > NGAS'

      do ilev=1,nlev
         dmf_h2o = vmr(ip_h2o,ilev)
         do kgas=1,ngas
            vmr(kgas,ilev) = vmr(kgas,ilev)/(1.+dmf_h2o)
         end do
      end do
      return
      end
