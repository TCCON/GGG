      subroutine legendre_poly_fit(mpts,npts,x,y,w,
     & nterm,jac,ip,wk,rnorm,a,b)
c  Fits a NTERM-term Legendre polynomial to the data [x,y] having weights w
c  Does this by solving the matrix equation  A.c=b
c  A(NPTS,NTERM) is the matrix of orthogonal Legendre polynomials
c  all scaled by the weights (W)
c  c(NTERM) is the vector of unknown coefficients
c  b(NPTS) is the measured values (y) scaled by the weights (W)
c
c  Inputs:
c                npts  I*4  number of data values
c             x(mpts)  R*4  vector of acscissae values
c             y(mpts)  R*4  vector of y-values
c             w(mpts)  R*4  vector of weights
c               nterm  I*4  Number of Legendre terms to be used
c     jac(mpts,nterm)  R*4  workspace array (Jacobian matrix)
c           ip(nterm)  I*4  workspace (for hfti)
c           wk(nterm)  R*4  workspace (for hfti)
c
c  Outputs:
c            y(nterm)  R*4  Legendre coefficient values
c               rnorm  R*4  Euclidean (L2) norm of the residual vector
c                 a b  R*4  Coefficents that transform x-coordinates 
c                           to zz=a*x+b that vary -1 to +1
c
c  Notes:
c        On output, the y-values are overwritten
c        The x-values have been normalized to fall into the range -1<x<+1
c
c  The first 6 Legendre polynomials are
c         L0(x) =  1
c         L1(x) =  x
c         L2(x) = (3x^2-1)/2
c         L3(x) = (5x^3-3x)/2
c         L4(x) = (35x^4-30x^2+3)/8
c         L5(x) = (63x^5-70x^3+15x)/8
c
c  The advantage of using Legendre polynomials rather
c  than regular polynomials (a + b.x + c.x^2 + d.x^3 + )
c  is that the terms are virtually orthogonal. They are
c  completely orthogonal if:
c   - the x-values are chosen to correspond to
c     the roots of the Legendre polynomials
c   - the function is over-sampled and the weights are uniform.
c
c  The advantages of orthogonality are:
c  - fewer terms are needed to achieve a given accuracy
c  - dropping terms has little impact on the lower terms
c  - lower sensitivity to rounding error.

      implicit none
      integer i,j,nterm,mpts,npts,ip(nterm),krank
      real*4 x(mpts),y(mpts),w(mpts),jac(mpts,nterm),
     &  wk(nterm),tau,rnorm,a,b,zz,wmax

      tau=1.0e-6  ! condition number for rank-deficiency
      if(nterm.le.1) then
        write(*,*)'Calling Legendre with NTERM=',nterm
        stop 'Waste of time calling Legendre with so few terms'
      endif
c
      wmax=0.0
      do i=1,npts
        if(w(i).gt.wmax) wmax=w(i)
      end do
      if(wmax.eq.0.0) stop 'wmax < 0'

c Use recurrence relationship to populate JAC matrix
      a=2.0/(npts-1)
      b=float(npts+1)/(npts-1)
      do i=1,npts
         zz=a*x(i)-b
         jac(i,1)=w(i)/wmax
         jac(i,2)=zz*w(i)/wmax
         do j=3,nterm
           jac(i,j)=((2*j-3)*zz*jac(i,j-1)-(j-2)*jac(i,j-2))/(j-1)
         end do
c         write(*,*)i,(jac(i,j),j=1,nterm)
         y(i)=y(i)*w(i)/wmax  ! multiply y-values by weights
      end do
c
c  Solve the matrix equation J.c=y by QR decomposition
      call shfti(jac,mpts,npts,nterm,y,npts,
     & 1,tau,krank,rnorm,wk,ip)
      if(krank.lt.nterm) then
         write(*,*)'krank,nterm=',krank,nterm
         stop 'Legendre polynomial matrix rank-deficient'
      endif

      return
      end

      subroutine legendre_poly_eval(nterm,a,b,c,npts,x,y)
c  Evaluates the Legendre polynomial defined by the coefficients C(nterm)
c  at the x-values provided.
c
c   y(i) = SUM_over_k  C(k)*L(k,x(i))
c
c  Inputs:
c            nterm  I*4  Number of Legendre terms to be used
c              a b  R*4  zz=a.x+b
c         c(nterm)  R*4  Legendre coefficients
c             npts  I*4  number of data values
c          x(mpts)  R*4  vector of abscissae values
c
c  Outputs:
c          y(mpts)  R*4  polynomial values

      implicit none
      integer i,j,nterm,npts
      real*4 x(npts),y(npts),c(nterm),yy,
     & bfunj,bfjm1,bfjm2,a,b,zz

c  Use recurrence relationship to evaluate Legendre basis functions (bfunj)
c  Then multiply each basis function by the appropriate coefficient.
      do i=1,npts
         zz=a*x(i)-b
         bfjm2=1.0
         bfjm1=zz
         yy=c(1)+c(2)*zz
         do j=3,nterm
            bfunj=((2*j-3)*zz*bfjm1-(j-2)*bfjm2)/(j-1)
            yy=yy+c(j)*bfunj
            bfjm2=bfjm1
            bfjm1=bfunj
         end do
         y(i)=yy
      end do
      return
      end
