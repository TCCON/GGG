      function  polyinterp(vin,nterm,xx)
c
c  Computes vin(xx) by polynomical interpolation
c  Special case of Legendre's equation for equally-spaced x-values
c  xx must be < 0.5 (i.e. vin(0) is the nearest point to the desired x-value).
c  Calling program must check vin(-nterm:nterm) doesn't exceed array bounds.
c  -ve NTERM forces the general code (useful for testing).
c  
      implicit none
      integer*4 nterm,k
      real*4 vin(-nterm:nterm),polyinterp,tt,cc,vv,xx
      real*4 fm4,fm3,fm2,fm1,fp1,fp2,fp3,fp4
    
c  For the most commonly used NTERM values (0-4), use the loop-unrolled
c  versions of the code (50% faster or 2/3 the time for NTERM=4)
      if(nterm.eq.0) then
         polyinterp=vin(0)
      elseif(nterm.eq.1) then  ! Quadratic Inperpolation
         fm1=xx-1
         fp1=xx+1
         polyinterp=fm1*fp1*
     &     (-vin(0) + xx*(
     &      +( vin(-1)/fp1 + vin(1)/fm1 )/2
     &     ))
      elseif(nterm.eq.2) then  ! Quartic Inperpolation
         fm2=xx-2
         fm1=xx-1
         fp1=xx+1
         fp2=xx+2
         polyinterp=fm2*fm1*fp1*fp2*
     &    ( +vin(0)/4 + xx*(
     &      (-( vin(-1)/fp1 + vin(1)/fm1 ))/6
     &      +( vin(-2)/fp2 + vin(2)/fm2 )/24
     &      ))
      elseif(nterm.eq.3) then  ! Sextic Inperpolation
         fm3=xx-3
         fm2=xx-2
         fm1=xx-1
         fp1=xx+1
         fp2=xx+2
         fp3=xx+3
         polyinterp=fm3*fm2*fm1*fp1*fp2*fp3*
     &     ((-vin(0))/36 + xx*(
     &      +( vin(-1)/fp1 + vin(1)/fm1 )/48
     &      -( vin(-2)/fp2 + vin(2)/fm2 )/120
     &      +( vin(-3)/fp3 + vin(3)/fm3 )/720
     &     ))
      elseif(nterm.eq.4) then  ! Octic Inperpolation
         fm4=xx-4
         fm3=xx-3
         fm2=xx-2
         fm1=xx-1
         fp1=xx+1
         fp2=xx+2
         fp3=xx+3
         fp4=xx+4
         polyinterp=fm4*fm3*fm2*fm1*fp1*fp2*fp3*fp4*
     &     (+vin(0)/576 + xx*(
     &      (-( vin(-1)/fp1 + vin(1)/fm1 ))/720
     &      +( vin(-2)/fp2 + vin(2)/fm2 )/1440
     &      -( vin(-3)/fp3 + vin(3)/fm3 )/5040
     &      +( vin(-4)/fp4 + vin(4)/fm4 )/40320
     &     ))
      else
c  The general code (below) works for any NTERM value but is slower.
         tt=1.0
         cc=1.0
         vv=0.0
         do k=iabs(nterm),1,-1
            cc=cc*(xx-k)*(xx+k)/(2*k-1)/(2*k)
            vv=vv+tt*(vin(-k)/(xx+k)+vin(k)/(xx-k))
            tt=(-tt)*(nterm+k)/(nterm-k+1)
         end do
         polyinterp=cc*(tt*vin(0)+xx*vv)
      endif

      return
      end

