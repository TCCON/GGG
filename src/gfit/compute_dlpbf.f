      subroutine compute_dlpbf(nmax,n,m,bf)
c  Computes Discrete Legendre Polynomial Basis Functions.
c  Based on the recursive method described in "Generation Scheme" 
c  section (p 758-759) of Neuman and Schonbach (1974).
c
c  BFs are normalized such that the first point has a value of +1.
c  Last point has a value of +1 or -1, depending whether it's even or odd.
c
c  Inputs:
c     nmax   I*4     Declared first dimension of array BF
c     n      I*4     Number of spectral points across window
c     m      I*4     Number of LP basis functions to compute
c
c  Output:
c   BF(nmax,m)  R*4     matrix of LP basis functions

c  k is the index over n (k=1,n).
c  i and j are indices over m (i=0,m)
c
c  Notes:
c
c  The Neuman and Schonbach implementation seems to assume that
c  the outer loop is over N (spectral point) and that the inner loop
c  is over M (order of the LP).
c
c  In my opinion, if N is much larger than M, then it is more efficient
c  to have the outer loop over M and the inner loop over N. This is
c  because there are more intermediate constants that depend on M
c  than depend on K=1,N. For this reason, and the fact that I want my
c  indices to run from 1,N (rather than 0,M), the implementation below
c  is slightly different from that in Neuman and Schonbach (1974).
c
c  Even though the first basis function (all 1's) may not
c  be used externally, it must still be computed since it
c  is used internally in the computation of the third BF.

      implicit none
      integer*4 nmax,n,m,i,j,k,nmo
      real*4 bf(nmax,m),cff,dd
      real*4 c0,c1,c2,d0,d2

      nmo=n-1
c  Compute first two LP basis functions
c  BF(*,1) is all ones,
      do k=1,n
         bf(k,1)=1.0
      end do

c  BF(*,2) is a straight line going from +1 to -1.
      if (m.ge.2) then
      do k=1,n
         bf(k,2)=dfloat(nmo-2*(k-1))/nmo
      end do
      endif

c  Now compute subsequent LP basis functions recursively.
c  So each basis function depends on the two previous ones
c  the coefficients c1, c2, c3, and BF(2).
      c0=0
      c1=nmo
      c2=nmo
      d0=nmo
      d2=nmo
      do i=3,m
         d0=d0+2
         d2=d2-2
         c0=c0+d0
         c1=c1+2*nmo
         c2=c2+d2
         do k=1,n
            bf(k,i) = ( c1*bf(k,2)*bf(k,i-1) - c0*bf(k,i-2) )/c2
         end do
      end do

c Renormalize DLPBFs such that their L2 norm is N
c Probably this re-normalization can be performed
c as part of the previous loops (future project).
c After re-normalization, the first basis function
c is still all ones, but the second (tilt) is reduced
c reduced by a factor 3(N-1)/(N+1), which for large N
c is approx SQRT(3).
      do i=1,m
         dd=sqrt(cff(n,i))
         do j=1,n
            bf(j,i) = bf(j,i)/dd 
         end do
      end do

      return
      end

      real*4 function cff(n,k)
c  Computes the L2 norm of an order-K, N-point, Discrete Legendre
c  Polynomial Basis Function (the dot product with itself).
c  Uses the first (un-numbered) equation in the Properties section
c  of Newmann and Schonbach [1974]. Adding one to all the indices,
c  because I like my loops to start at 1, this equation becomes
c       (N+K-1)! (N-K)! / ((N-1)!)^2 / (2*K-1)
c  This subroutine avoids the overflow associated with
c  taking actual factorials for large N. Instead, it
c  computes the ratios of the numerator and denominator
c  factorial terms together.
      integer n,k,j
c      cff=float(n)/(2*k-1)   ! L2-norm=1
      cff=1.0/(2*k-1)         ! L2-norm=N
      do j=1+n,k+n-1
         cff=cff*float(j)/(j-k)
      end do
      return
      end
