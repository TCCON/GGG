      subroutine symmetry(rbuf,symmi,symmh,np)
c
c  Computes the symmetry of an interferogram about the point RBUF(0)
c  and also about RBUF(0.5), the mid-point between RBUF(0) and RBUF(1).
c
c  INPUTS:
c     NP               I*4   # points over which symmetry is to be evaluated
c     RBUF(-np/2:np/2) R*4   section of raw interferogram
c
c  OUTPUTS:
c     SYMMI            R*4   interferogram symmetry about RBUF(0.0)
c     SYMMH            R*4   interferogram symmetry about RBUF(0.5)
c
c  Symmetry is evaluated by comparing the sums and differences of
c  points that are symmetrical about RBUF(0)
c      SYMMI = (SPI-SMI)/(SPI+SMI)
c
c  where SMI = SUM |RBUF(-i) - RBUF(+i)|
c        SPI = SUM |RBUF(-i) + RBUF(+i)|
c
c  If igram is perfectly symmetrical about RBUF(0),      SMI=0, SYMMI=+1
c  If igram is perfectly anti-symmetrical about RBUF(0), SPI=0, SYMMI=-1
c
c  Similary, the symmetry can also be defined about the half-point RBUF(0.5)
c  using  SMH = SUM |RBUF(-i+1) - RBUF(+i)|
c         SPH = SUM |RBUF(-i+1) + RBUF(+i)|
c       SYMMH = (SPH-SMH)/(SPH+SMH)

      implicit none

      integer*4 np,j
      real*4 rbuf(-(np/2):+np/2),tiny,
     & spi,smi,sph,smh,symmi,symmh,q,x,ww,pi
      parameter (tiny=1.E-36, pi=real(4.d0*datan(1.d0)))

      spi=0.0
      smi=0.0
      sph=0.0
      smh=0.0
      q=pi/float(np)
      do j=1,np/2
         x=q*float(j)
         ww=(5.0*cos(x)+cos(3.0*x))/6.0
         spi=spi+ww*abs(rbuf(-j)+rbuf(j))
         smi=smi+ww*abs(rbuf(-j)-rbuf(j))
         sph=sph+ww*abs(rbuf(-j+1)+rbuf(j))
         smh=smh+ww*abs(rbuf(-j+1)-rbuf(j))
      enddo
      symmi=(spi-smi)/(spi+smi+tiny)
      symmh=(sph-smh)/(sph+smh+tiny)
      return
      end
