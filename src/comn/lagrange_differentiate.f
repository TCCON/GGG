      subroutine lagrange_differentiate(np,vin,vout)
c  Computes the gradients of the function VIN.
c  Gradients are evaluated at the at the points themselves.
c
c  Inputs:
c            NP   I*4  Number of points in in/output vectors
c        VIN(NP)  R*4  Input vector
c
c  Outputs:
c       VOUT(NP)  R*4  Output vector (containing gradients)
c
c  Does this by fitting a Lagrange polynomial to the input vector,
c  differentiating analytically, then evaluating the differentials
c  at the equally-spaced abscissae.
c  Uses 3-point Lagrange interpolation for the edge points
c  Uses 4-point Lagrange interpolation for the penultimate points
c  Uses 5-point Lagrange interpolation for the rest of the window
c  The 5-point operator only has 4 terms, because the coefficient
c  of the VIN(i) term is zero. 
c
c  Gradients are in units of VIN per point spacing. So if you halve
c  the grid, the gradients will also halve.

      integer*4 np,i
      real*4 vin(np),vout(np)

      if(np.ge.4) then
         vout(1)=(-3*vin(1)+4*vin(2)-vin(3))/2                   ! 3-point
         vout(2)=(-2*vin(1)-3*vin(2)+6*vin(3)-vin(4))/6          ! 4-point
         do i=1+2,np-2
            vout(i)=(vin(i-2)-8*vin(i-1)+8*vin(i+1)-vin(i+2))/12 ! 5-point
         end do
         vout(np-1)=(2*vin(np)+3*vin(np-1)-6*vin(np-2)+vin(np-3))/6 ! 4-point
         vout(np)=(3*vin(np)-4*vin(np-1)+vin(np-2))/2            ! 3-point
      elseif(np.eq.3) then
         vout(1)=(-3*vin(1)+4*vin(2)-vin(3))/2
         vout(2)=(-vin(1)+vin(3))/2
         vout(3)=(3*vin(1)-4*vin(2)+vin(3))/2
      elseif(np.eq.2) then
         vout(1)=(-vin(1)+vin(2))/2
         vout(2)=(-vin(1)+vin(2))/2
      elseif(np.eq.1) then
         vout(1)=0.0
      endif
      return
      end

