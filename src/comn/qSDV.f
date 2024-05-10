        subroutine qSDV(sg0,cte,Gam0,Gam2,Shift0,Shift2,
     &sg,Y_LM,N,g,K)
!-------------------------------------------------
!       "qSDV": quadratic-Speed-Dependent Voigt
!       Subroutine to Compute the complex normalized spectral shape of an 
!       isolated line by the qSDV model
!
!       Input/Output Parameters of Routine (Arguments or Common)
!       ---------------------------------
!       sg0     : Unperturbed line position in cm-1 (Input).
!       cte     : Doppler HWHM in cm-1 (Input)
!       Gam0    : Speed-averaged line-width in cm-1 (Input). 
!       Gam2    : Speed dependence of the line-width in cm-1 (Input).
!       Shift0  : Speed-averaged line-shift in cm-1 (Input).
!       Shift2  : Speed dependence of the line-shift in cm-1 (Input)
!       sg      : Current WaveNumber of the Computation in cm-1 (Input).
!       N       : Number of grid points
!       g       : grid 
!
!       Output Quantities (through Common Statements)
!       -----------------
!       LS_qSDV_R: Real part of the normalized spectral shape (cm)
!       LS_qSDV_I: Imaginary part of the normalized spectral shape (cm)
!
!       Called Routines: 'CPF1' (Complex Probability Function)
!       ---------------  'CPF3' (Complex Probability Function for the region 3)
!
!       Called By: Main Program
!       ---------
!
!     Double Precision Version
!
!-------------------------------------------------
        implicit none
        double precision sg0
        double precision Gam0,Gam2,Shift0,Shift2
        double precision sg,g
        double precision pi,rpi
        REAL*4 cte
        INTEGER*4 N,I
        double precision xz1,xz2,yz1,yz2,xXb,yXb
        double precision wr1,wi1,wr2,wi2,wrb,wib
        double precision SZ1,SZ2,DSZ,SZmx,SZmn
        double precision LS_qSDV_R,LS_qSDV_I,Y_LM!,K(0:N-1)
        Real*4 K(0:N-1)
        double complex c0,c2,c0t,c2t
        double complex X,Y,iz,Z1,Z2
        double complex Aterm,LS_qSDV
!
!-------------------------------------------------
!
        cte=1/cte
        pi=3.141592654d0
        rpi=dsqrt(pi)
        iz=dcmplx(0.d0,1.d0)
! Calculating the different parameters 
        c0=dcmplx(Gam0,-Shift0)
        c2=dcmplx(Gam2,-Shift2)
        c0t=(c0-1.5d0*c2)
        c2t=c2

C Start loop over all grid points        
        DO I=0,N-1
        if (cdabs(C2t).eq.0.d0) go to 110
        Y=1.d0/((2.d0*cte*C2t))**2
        X=(iz*(sg-sg0)+c0t)/c2t

        if (cdabs(X).le.3.d-8*cdabs(Y)) go to 120
        if (cdabs(Y).le.1.d-15*cdabs(X)) go to 140
! calculating Z1 and Z2
        Z1=cdsqrt(X+Y)-cdsqrt(Y)
        Z2=Z1+2.d0*cdsqrt(Y)
! calculating the real and imaginary parts of Z1 and Z2
        xZ1=-dimag(Z1)
        yZ1=dreal(Z1)
        xZ2=-dimag(Z2)
        yZ2=dreal(Z2)
! check if Z1 and Z2 are close to each other
        SZ1=dsqrt(xZ1*xZ1+yZ1*yZ1)
        SZ2=dsqrt(xZ2*xZ2+yZ2*yZ2)
        DSZ=dabs(SZ1-SZ2)
        SZmx=dmax1(SZ1,SZ2)
        SZmn=dmin1(SZ1,SZ2)
! when Z1 and Z2 are close to each other, ensure that they are in 
! the same interval of CPF1 
        if (DSZ.le.1.d0.and.SZmx.gt.8.d0.and.SZmn.le.8.d0) then
        Call CPF3_3(xZ1,yZ1,WR1,WI1) 
        Call CPF3_3(xZ2,yZ2,WR2,WI2) 
        else
        Call CPF1(xZ1,yZ1,WR1,WI1) 
        Call CPF1(xZ2,yZ2,WR2,WI2) 
        endif
! calculating the A term of the profile
!        Aterm=rpi*cte*(dcmplx(wr1,wi1)-dcmplx(wr2,wi2))
        Aterm=dcmplx(WR1,WI1)-dcmplx(WR2,WI2)
        go to 10
! when C2t=0 normal Voigt calculation
110   continue
        Z1=(iz*(sg-sg0)+C0t)*cte
        xZ1=-dimag(Z1)
        yZ1=dreal(Z1)
        Call CPF1(xZ1,yZ1,WR1,WI1)
!        Aterm=rpi*cte*dcmplx(WR1,WI1)
        Aterm=dcmplx(WR1,WI1) 
       go to 10
! when abs(Y) is much larger than abs(X)
120   continue
        Z1=(iz*(sg-sg0)+C0t)*cte
        Z2=cdsqrt(X+Y)+cdsqrt(Y)
        xZ1=-dimag(z1)
        yZ1=dreal(z1)
        xZ2=-dimag(z2)
        yZ2=dreal(z2)
        Call CPF1(xZ1,yZ1,WR1,WI1)
        Call CPF1(xZ2,yZ2,WR2,WI2) 
!        Aterm=rpi*cte*(dcmplx(WR1,WI1)-dcmplx(WR2,WI2))
        Aterm=dcmplx(WR1,WI1)-dcmplx(WR2,WI2)
        go to 10
! when abs(X) is much larger than abs(Y)
140   continue
        if (cdabs(cdsqrt(X)).le.4.d3) then
          xXb=-dimag(cdsqrt(X))
          yXb=dreal(cdsqrt(X))
          Call CPF1(xXb,yXb,WRb,WIb) 
          Aterm=(2.d0*rpi/C2t)*(1.d0/rpi-cdsqrt(X)*dcmplx(WRb,WIb))
          write(*,*)'X>>Y',sg0,sg
! and when abs(X) is much larger than 1
        else
          Aterm=(1.d0/C2t)*(1.d0/X-1.5d0/(X**2))
        endif
!
10    continue
!
!       LS_qSDV=(1.d0/pi)*Aterm
       LS_qSDV=Aterm
       LS_qSDV_R=dreal(LS_qSDV)
       LS_qSDV_I=dimag(LS_qSDV)
       K(I)=(LS_qSDV_R-(Y_LM*LS_qSDV_I))
!       K(I)=LS_qSDV_R
       sg=sg+g
      ENDDO   
!      write(*,*)'From SD:',K,'Y=',Y_LM,'R=',LS_qSDV_R,'I=',LS_qSDV_I
      Return
      End Subroutine qSDV

            Subroutine CPF1(X,Y,WR,WI)
!-------------------------------------------------
! "CPF": Complex Probability Function
! .........................................................
!         .       Subroutine to Compute the Complex       .
!         .        Probability Function W(z=X+iY)         .
!         .     W(z)=exp(-z**2)*Erfc(-i*z) with Y>=0      .
!         .    Which Appears when Convoluting a Complex   .
!         .     Lorentzian Profile by a Gaussian Shape    .
!         .................................................
!
!             WR : Real Part of W(z)
!             WI : Imaginary Part of W(z)
!
! This Routine was Taken from the Paper by J. Humlicek, which 
! is Available in Page 309 of Volume 21 of the 1979 Issue of
! the Journal of Quantitative Spectroscopy and Radiative Transfer
! Please Refer to this Paper for More Information
!
! Accessed Files:  None
! --------------
!
! Called Routines: None                               
! ---------------                                 
!
! Called By: 'CompAbs' (COMPute ABSorpton)
! ---------
!
! Double Precision Version
!
!-------------------------------------------------
!      
      Implicit None
      Integer I
      double complex zm1,zm2,zterm,zsum,zone,zi
      Double Precision X,Y,WR,WI
      Double Precision T,U,S,Y1,Y2,Y3,R,R2,D,D1,D2,D3,D4
      Double Precision TT(15),pipwoeronehalf
!      
      Dimension T(6),U(6),S(6)
      Data T/.314240376d0,.947788391d0,1.59768264d0,2.27950708d0
     ,        ,3.02063703d0,3.8897249d0/
      Data U/1.01172805d0,-.75197147d0,1.2557727d-2,1.00220082d-2
     ,        ,-2.42068135d-4,5.00848061d-7/
      Data S/1.393237d0,.231152406d0,-.155351466d0,6.21836624d-3
     ,        ,9.19082986d-5,-6.27525958d-7/
      Data zone,zi/(1.d0,0.D0),(0.d0,1.D0)/
      data tt/0.5d0,1.5d0,2.5d0,3.5d0,4.5d0,5.5d0,6.5d0,7.5d0,8.5d0,
     ,        9.5d0,10.5d0,11.5d0,12.5d0,13.5d0,14.5d0/
      data pipwoeronehalf/0.564189583547756d0/

! new Region 3
        if(dsqrt(x*x+y*Y).gt.8.D0)then
        zm1=zone/dcmplx(x,y)
        zm2=zm1*zm1
        zsum=zone
        zterm=zone
        do i=1,15
        zterm=zterm*zm2*tt(i)
        zsum=zsum+zterm
        end do
        zsum=zsum*zi*zm1*pipwoeronehalf
        wr=dreal(zsum)
        wi=dimag(zsum)
        return
        end if
!
      WR=0.d0
      WI=0.d0
      Y1=Y+1.5d0
      Y2=Y1*Y1
      If( (Y.GT.0.85d0) .OR. (DABS(X).LT.(18.1d0*Y+1.65d0)) )GoTo 2
!
!       Region 2
!
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
!
!       Region 1
!
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
      End Subroutine CPF1

            Subroutine CPF3_3(X,Y,WR,WI)
!-------------------------------------------------
! "CPF": Complex Probability Function
! .........................................................
!         .       Subroutine to Compute the Complex       .
!         .        Probability Function W(z=X+iY)         .
!         .     W(z)=exp(-z**2)*Erfc(-i*z) with Y>=0      .
!         .    Which Appears when Convoluting a Complex   .
!         .     Lorentzian Profile by a Gaussian Shape    .
!         .................................................
!
!             WR : Real Part of W(z)
!             WI : Imaginary Part of W(z)
!
! This Routine takes into account the region 3 only, i.e. when sqrt(x**2+y**2)>8. 
!
! Accessed Files:  None
! --------------
!
! Called Routines: None                               
! ---------------                                 
!
! Called By: 'pCqSDHC'
! ---------
!
! Double Precision Version
! 
!-------------------------------------------------
!      
      Implicit None
      Integer I
      double complex zm1,zm2,zterm,zsum,zone,zi
      Double Precision X,Y,WR,WI
      Double Precision TT(15),pipwoeronehalf
!      
        Data zone,zi/(1.d0,0.D0),(0.d0,1.D0)/
        data tt/0.5d0,1.5d0,2.5d0,3.5d0,4.5d0,5.5d0,6.5d0,7.5d0,8.5d0,
     ,        9.5d0,10.5d0,11.5d0,12.5d0,13.5d0,14.5d0/
        data pipwoeronehalf/0.564189583547756d0/

! Region 3
        zm1=zone/dcmplx(x,y)
        zm2=zm1*zm1
        zsum=zone
        zterm=zone
        do i=1,15
        zterm=zterm*zm2*tt(i)
        zsum=zsum+zterm
        end do
        zsum=zsum*zi*zm1*pipwoeronehalf
        wr=dreal(zsum)
        wi=dimag(zsum)
       return
      End Subroutine CPF3_3

