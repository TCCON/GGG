C*********************************************************************
      Subroutine CPF(X,Y,WR,WI)
C*********************************************************************
C	"CPF": Complex Probability Function
C	.................................................
C	.	Subroutine to Compute the Complex           .
C	.	Probability Function W(z=X+iY)              .
C	.	W(z)=exp(-z**2)*Erfc(-i*z) with Y>=0        .
C	.	Which Appears when Convoluting a Complex    .
C	.	Lorentzian Profile by a Gaussian Shape      .
C	.................................................
C
C	WR	 : Real Part of W(z)
C	WI	 : Imaginary Part of W(z)
C
C	This Routine was Taken from the Paper by J. Humlicek, which 
C	is Available in Page 309 of Volume 21 of the 1979 Issue of
C	the Journal of Quantitative Spectroscopy and Radiative Transfer
C	Please Refer to this Paper for More Information
C
C	Accessed Files:  None
C	--------------
C
C	Called Routines: None                               
C	---------------                                 
C
C	Called By: 'CompAbs' (COMPute ABSorpton)
C	---------
C
C	Double Precision Version
C
C	J.-M. Hartmann, last change 06 March 1997
C*********************************************************************
C      
      Implicit None
      Integer I
      Double Precision X,Y,WR,WI
      Double Precision T,U,S,Y1,Y2,Y3,R,R2,D,D1,D2,D3,D4
C      
      Dimension T(6),U(6),S(6)
      Data T/.314240376d0,.947788391d0,1.59768264d0,2.27950708d0
     ,        ,3.02063703d0,3.8897249d0/
      Data U/1.01172805d0,-.75197147d0,1.2557727d-2,1.00220082d-2
     ,        ,-2.42068135d-4,5.00848061d-7/
      Data S/1.393237d0,.231152406d0,-.155351466d0,6.21836624d-3
     ,        ,9.19082986d-5,-6.27525958d-7/
      WR=0.d0
      WI=0.d0
      Y1=Y+1.5d0
      Y2=Y1*Y1
      If( (Y.GT.0.85d0) .OR. (DABS(X).LT.(18.1d0*Y+1.65d0)) )GoTo 2
C
C       Region II
C
      If( DABS(X).LT.12.d0 )WR=DEXP(-X*X)
      Y3=Y+3.d0
      Do 1 I=1,6
      R=X-T(I)
      R2=R*R
      D=1.d0/(R2+Y2)
      D1=Y1*D
      D2=R*D
      WR=WR+Y*(U(I)*(R*D2-1.5d0*D1)+S(I)*Y3*D2)/(R2+2.25d0)
      R=X+T(I)
      R2=R*R
      D=1.d0/(R2+Y2)
      D3=Y1*D
      D4=R*D
      WR=WR+Y*(U(I)*(R*D4-1.5d0*D3)-S(I)*Y3*D4)/(R2+2.25d0)
      WI=WI+U(I)*(D2+D4)+S(I)*(D1-D3)
 1    Continue  
      Return
C
C       Region I
C
 2    Continue
      Do 3 I=1,6
      R=X-T(I)
      D=1.d0/(R*R+Y2)
      D1=Y1*D
      D2=R*D
      R=X+T(I)
      D=1.d0/(R*R+Y2)
      D3=Y1*D
      D4=R*D
      WR=WR+U(I)*(D1+D3)-S(I)*(D2-D4)
      WI=WI+U(I)*(D2+D4)+S(I)*(D1-D3)
 3    Continue  
      Return
      End Subroutine CPF
