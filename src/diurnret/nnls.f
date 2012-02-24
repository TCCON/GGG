C     SUBROUTINE NNLS  (A,MDA,M,N,B,X,RNORM,W,ZZ,INDEX,MODE)
C     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUNE 15
C     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974
C   
C         **********   NONNEGATIVE LEAST SQUARES   **********   
C   
C     GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE AN
C     N-VECTOR, X, WHICH SOLVES THE LEAST SQUARES PROBLEM   
C   
C                      A * X = B  SUBJECT TO X .GE. 0   
C   
C     A(),MDA,M,N     MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE   
C                     ARRAY, A().   ON ENTRY A() CONTAINS THE M BY N    
C                     MATRIX, A.           ON EXIT A() CONTAINS 
C                     THE PRODUCT MATRIX, Q*A , WHERE Q IS AN   
C                     M BY M ORTHOGONAL MATRIX GENERATED IMPLICITLY BY  
C                     THIS SUBROUTINE.  
C     B()     ON ENTRY B() CONTAINS THE M-VECTOR, B.   ON EXIT B() CON- 
C             TAINS Q*B.
C     X()     ON ENTRY X() NEED NOT BE INITIALIZED.  ON EXIT X() WILL   
C             CONTAIN THE SOLUTION VECTOR.  
C     RNORM   ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE  
C             RESIDUAL VECTOR.  
C     W()     AN N-ARRAY OF WORKING SPACE.  ON EXIT W() WILL CONTAIN    
C             THE DUAL SOLUTION VECTOR.   W WILL SATISFY W(I) = 0.  
C             FOR ALL I IN SET P  AND W(I) .LE. 0. FOR ALL I IN SET Z   
C     ZZ()     AN M-ARRAY OF WORKING SPACE.     
C     INDEX()     AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N.
C                 ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS    
C                 P AND Z AS FOLLOWS..  
C   
C                 INDEX(1)   THRU INDEX(NSETP) = SET P.     
C                 INDEX(IZ1) THRU INDEX(IZ2)   = SET Z.     
C                 IZ1 = NSETP + 1 = NPP1
C                 IZ2 = N   
C     MODE    THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING 
C             MEANINGS. 
C             1     THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY.
C             2     THE DIMENSIONS OF THE PROBLEM ARE BAD.  
C                   EITHER M .LE. 0 OR N .LE. 0.
C             3    ITERATION COUNT EXCEEDED.  MORE THAN 3*N ITERATIONS. 
C   
      SUBROUTINE NNLS (A,MDA,M,N,B,X,RNORM,W,ZZ,INDEX,MODE) 
      implicit none
      integer mda,m,n,mode,iter,i,itmax,iz1,iz2,nsetp,npp1,iz,j,l,izmax,
     & jz,jj,next,ii,ip
      real rnorm,zero,one,two,factor,ztest,alpha,t,cc,ss,diff,sm,wmax,
     & asave,up,dummy,unorm
      REAL A(MDA,N), B(M), X(N), W(N), ZZ(M)   
      INTEGER INDEX(N)  
      ZERO=0.   
      ONE=1.
      TWO=2.
      FACTOR=0.01   
C   
      MODE=1
      IF (M.GT.0.AND.N.GT.0) GO TO 10   
      MODE=2
      RETURN
   10 ITER=0
      ITMAX=3*N 
C   
C                    INITIALIZE THE ARRAYS INDEX() AND X(). 
C   
          DO 20 I=1,N   
          X(I)=ZERO     
   20     INDEX(I)=I    
C   
      IZ2=N 
      IZ1=1 
      NSETP=0   
      NPP1=1
C                             ******  MAIN LOOP BEGINS HERE  ******     
   30 CONTINUE  
C                  QUIT IF ALL COEFFICIENTS ARE ALREADY IN THE SOLUTION.
C                        OR IF M COLS OF A HAVE BEEN TRIANGULARIZED.    
C   
      IF (IZ1.GT.IZ2.OR.NSETP.GE.M) GO TO 350   
C   
C         COMPUTE COMPONENTS OF THE DUAL (NEGATIVE GRADIENT) VECTOR W().
C   
          DO 50 IZ=IZ1,IZ2  
          J=INDEX(IZ)   
          SM=ZERO   
              DO 40 L=NPP1,M
   40         SM=SM+A(L,J)*B(L)     
   50     W(J)=SM   
C                                   FIND LARGEST POSITIVE W(J). 
   60 WMAX=ZERO 
          DO 70 IZ=IZ1,IZ2  
          J=INDEX(IZ)   
          IF (W(J).LE.WMAX) GO TO 70
          WMAX=W(J)     
          IZMAX=IZ  
   70     CONTINUE  
C   
C             IF WMAX .LE. 0. GO TO TERMINATION.
C             THIS INDICATES SATISFACTION OF THE KUHN-TUCKER CONDITIONS.
C   
      IF (WMAX) 350,350,80  
   80 IZ=IZMAX  
      J=INDEX(IZ)   
C   
C     THE SIGN OF W(J) IS OK FOR J TO BE MOVED TO SET P.    
C     BEGIN THE TRANSFORMATION AND CHECK NEW DIAGONAL ELEMENT TO AVOID  
C     NEAR LINEAR DEPENDENCE.   
C   
      ASAVE=A(NPP1,J)   
      CALL H12 (1,NPP1,NPP1+1,M,A(1,J),1,UP,DUMMY,1,1,0)    
      UNORM=ZERO
      IF (NSETP.EQ.0) GO TO 100     
          DO 90 L=1,NSETP   
   90     UNORM=UNORM+A(L,J)**2     
  100 UNORM=SQRT(UNORM) 
      IF (DIFF(UNORM+ABS(A(NPP1,J))*FACTOR,UNORM)) 130,130,110  
C   
C     COL J IS SUFFICIENTLY INDEPENDENT.  COPY B INTO ZZ, UPDATE ZZ AND 
C     SOLVE FOR ZTEST ( = PROPOSED NEW VALUE FOR X(J) ).    
C   
  110     DO 120 L=1,M  
  120     ZZ(L)=B(L)    
      CALL H12 (2,NPP1,NPP1+1,M,A(1,J),1,UP,ZZ,1,1,1)   
      ZTEST=ZZ(NPP1)/A(NPP1,J)  
C   
C                                     SEE IF ZTEST IS POSITIVE  
C   
      IF (ZTEST) 130,130,140
C   
C     REJECT J AS A CANDIDATE TO BE MOVED FROM SET Z TO SET P.  
C     RESTORE A(NPP1,J), SET W(J)=0., AND LOOP BACK TO TEST DUAL
C     COEFFS AGAIN.     
C   
  130 A(NPP1,J)=ASAVE   
      W(J)=ZERO 
      GO TO 60  
C   
C     THE INDEX  J=INDEX(IZ)  HAS BEEN SELECTED TO BE MOVED FROM
C     SET Z TO SET P.    UPDATE B,  UPDATE INDICES,  APPLY HOUSEHOLDER  
C     TRANSFORMATIONS TO COLS IN NEW SET Z,  ZERO SUBDIAGONAL ELTS IN   
C     COL J,  SET W(J)=0.   
C   
  140     DO 150 L=1,M  
  150     B(L)=ZZ(L)    
C   
      INDEX(IZ)=INDEX(IZ1)  
      INDEX(IZ1)=J  
      IZ1=IZ1+1 
      NSETP=NPP1
      NPP1=NPP1+1   
C   
      IF (IZ1.GT.IZ2) GO TO 170     
          DO 160 JZ=IZ1,IZ2 
          JJ=INDEX(JZ)  
  160     CALL H12 (2,NSETP,NPP1,M,A(1,J),1,UP,A(1,JJ),1,MDA,1) 
  170 CONTINUE  
C   
      IF (NSETP.EQ.M) GO TO 190     
          DO 180 L=NPP1,M   
  180     A(L,J)=ZERO   
  190 CONTINUE  
C   
      W(J)=ZERO 
C                                SOLVE THE TRIANGULAR SYSTEM.   
C                                STORE THE SOLUTION TEMPORARILY IN ZZ().
      ASSIGN 200 TO NEXT
      GO TO 400 
  200 CONTINUE  
C   
C                       ******  SECONDARY LOOP BEGINS HERE ******   
C   
C                          ITERATION COUNTER.   
C   
  210 ITER=ITER+1   
      IF (ITER.LE.ITMAX) GO TO 220  
      MODE=3
c      WRITE (6,440) iter     
      GO TO 350 
  220 CONTINUE  
C   
C                    SEE IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE.    
C                                  IF NOT COMPUTE ALPHA.    
C   
      ALPHA=TWO 
          DO 240 IP=1,NSETP 
          L=INDEX(IP)   
          IF (ZZ(IP)) 230,230,240   
C   
  230     T=-X(L)/(ZZ(IP)-X(L))     
          IF (ALPHA.LE.T) GO TO 240 
          ALPHA=T   
          JJ=IP 
  240     CONTINUE  
C   
C          IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE THEN ALPHA WILL   
C          STILL = 2.    IF SO EXIT FROM SECONDARY LOOP TO MAIN LOOP.   
C   
      IF (ALPHA.EQ.TWO) GO TO 330   
C   
C          OTHERWISE USE ALPHA WHICH WILL BE BETWEEN 0. AND 1. TO   
C          INTERPOLATE BETWEEN THE OLD X AND THE NEW ZZ.    
C   
          DO 250 IP=1,NSETP 
          L=INDEX(IP)   
  250     X(L)=X(L)+ALPHA*(ZZ(IP)-X(L)) 
C   
C        MODIFY A AND B AND THE INDEX ARRAYS TO MOVE COEFFICIENT I  
C        FROM SET P TO SET Z.   
C   
      I=INDEX(JJ)   
  260 X(I)=ZERO 
C   
      IF (JJ.EQ.NSETP) GO TO 290    
      JJ=JJ+1   
          DO 280 J=JJ,NSETP 
          II=INDEX(J)   
          INDEX(J-1)=II 
          CALL G1 (A(J-1,II),A(J,II),CC,SS,A(J-1,II))   
          A(J,II)=ZERO  
              DO 270 L=1,N  
              IF (L.NE.II) CALL G2 (CC,SS,A(J-1,L),A(J,L))  
  270         CONTINUE  
  280     CALL G2 (CC,SS,B(J-1),B(J))   
  290 NPP1=NSETP
      NSETP=NSETP-1     
      IZ1=IZ1-1 
      INDEX(IZ1)=I  
C   
C        SEE IF THE REMAINING COEFFS IN SET P ARE FEASIBLE.  THEY SHOULD
C        BE BECAUSE OF THE WAY ALPHA WAS DETERMINED.
C        IF ANY ARE INFEASIBLE IT IS DUE TO ROUND-OFF ERROR.  ANY   
C        THAT ARE NONPOSITIVE WILL BE SET TO ZERO   
C        AND MOVED FROM SET P TO SET Z. 
C   
          DO 300 JJ=1,NSETP 
          I=INDEX(JJ)   
          IF (X(I)) 260,260,300     
  300     CONTINUE  
C   
C         COPY B( ) INTO ZZ( ).  THEN SOLVE AGAIN AND LOOP BACK.
C   
          DO 310 I=1,M  
  310     ZZ(I)=B(I)    
      ASSIGN 320 TO NEXT
      GO TO 400 
  320 CONTINUE  
      GO TO 210 
C                      ******  END OF SECONDARY LOOP  ******
C   
  330     DO 340 IP=1,NSETP 
          I=INDEX(IP)   
  340     X(I)=ZZ(IP)   
C        ALL NEW COEFFS ARE POSITIVE.  LOOP BACK TO BEGINNING.  
      GO TO 30  
C   
C                        ******  END OF MAIN LOOP  ******   
C   
C                        COME TO HERE FOR TERMINATION.  
C                     COMPUTE THE NORM OF THE FINAL RESIDUAL VECTOR.    
C   
  350 SM=ZERO   
      IF (NPP1.GT.M) GO TO 370  
          DO 360 I=NPP1,M   
  360     SM=SM+B(I)**2 
      GO TO 390 
  370     DO 380 J=1,N  
  380     W(J)=ZERO     
  390 RNORM=SQRT(SM)    
      RETURN
C   
C     THE FOLLOWING BLOCK OF CODE IS USED AS AN INTERNAL SUBROUTINE     
C     TO SOLVE THE TRIANGULAR SYSTEM, PUTTING THE SOLUTION IN ZZ().     
C   
  400     DO 430 L=1,NSETP  
          IP=NSETP+1-L  
          IF (L.EQ.1) GO TO 420     
              DO 410 II=1,IP
  410         ZZ(II)=ZZ(II)-A(II,JJ)*ZZ(IP+1)   
  420     JJ=INDEX(IP)  
  430     ZZ(IP)=ZZ(IP)/A(IP,JJ)    
      GO TO NEXT, (200,320) 
  440 FORMAT ('NNLS QUITTING AFTER',i4,' iterations')   
      END   


C     SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)  
C     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUN 12 
C     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974
C   
C     CONSTRUCTION AND/OR APPLICATION OF A SINGLE   
C     HOUSEHOLDER TRANSFORMATION..     Q = I + U*(U**T)/B   
C   
C     MODE    = 1 OR 2   TO SELECT ALGORITHM  H1  OR  H2 .  
C     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT. 
C     L1,M   IF L1 .LE. M   THE TRANSFORMATION WILL BE CONSTRUCTED TO   
C            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.   IF L1 GT. M     
C            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
C     U(),IUE,UP    ON ENTRY TO H1 U() CONTAINS THE PIVOT VECTOR.   
C                   IUE IS THE STORAGE INCREMENT BETWEEN ELEMENTS.  
C                                       ON EXIT FROM H1 U() AND UP  
C                   CONTAIN QUANTITIES DEFINING THE VECTOR U OF THE     
C                   HOUSEHOLDER TRANSFORMATION.   ON ENTRY TO H2 U()    
C                   AND UP SHOULD CONTAIN QUANTITIES PREVIOUSLY COMPUTED
C                   BY H1.  THESE WILL NOT BE MODIFIED BY H2.   
C     C()    ON ENTRY TO H1 OR H2 C() CONTAINS A MATRIX WHICH WILL BE   
C            REGARDED AS A SET OF VECTORS TO WHICH THE HOUSEHOLDER  
C            TRANSFORMATION IS TO BE APPLIED.  ON EXIT C() CONTAINS THE 
C            SET OF TRANSFORMED VECTORS.
C     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C().  
C     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C().  
C     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. IF NCV .LE. 0  
C            NO OPERATIONS WILL BE DONE ON C(). 
C   
      SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)  
      implicit none
      integer mode,lpivot,l1,m,iue,ice,icv,ncv,j,i,i2,i3,i4,incr
      real up,one,cl,clinv,sm1
      REAL U(IUE,M), C(*)  
      DOUBLE PRECISION SM,B 
      ONE=1.
C   
      IF (0.GE.LPIVOT.OR.LPIVOT.GE.L1.OR.L1.GT.M) RETURN    
      CL=ABS(U(1,LPIVOT))   
      IF (MODE.EQ.2) GO TO 60   
C                            ****** CONSTRUCT THE TRANSFORMATION. ******
          DO 10 J=L1,M  
   10     CL=AMAX1(ABS(U(1,J)),CL)  
      IF (CL) 130,130,20
   20 CLINV=ONE/CL  
      SM=(DBLE(U(1,LPIVOT))*CLINV)**2   
          DO 30 J=L1,M  
   30     SM=SM+(DBLE(U(1,J))*CLINV)**2 
C                              CONVERT DBLE. PREC. SM TO SNGL. PREC. SM1
      SM1=SM
      CL=CL*SQRT(SM1)   
      IF (U(1,LPIVOT)) 50,50,40     
   40 CL=-CL
   50 UP=U(1,LPIVOT)-CL 
      U(1,LPIVOT)=CL    
      GO TO 70  
C            ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******
C   
   60 IF (CL) 130,130,70
   70 IF (NCV.LE.0) RETURN  
      B=DBLE(UP)*U(1,LPIVOT)
C                       B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.
C   
      IF (B) 80,130,130 
   80 B=ONE/B   
      I2=1-ICV+ICE*(LPIVOT-1)   
      INCR=ICE*(L1-LPIVOT)  
          DO 120 J=1,NCV
          I2=I2+ICV     
          I3=I2+INCR    
          I4=I3 
          SM=C(I2)*DBLE(UP) 
              DO 90 I=L1,M  
              SM=SM+C(I3)*DBLE(U(1,I))  
   90         I3=I3+ICE 
          IF (SM) 100,120,100   
  100     SM=SM*B   
          C(I2)=C(I2)+SM*DBLE(UP)   
              DO 110 I=L1,M 
              C(I4)=C(I4)+SM*DBLE(U(1,I))   
  110         I4=I4+ICE 
  120     CONTINUE  
  130 RETURN
      END   


      SUBROUTINE G1 (A,B,COS,SIN,SIG)   
C     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUN 12 
C     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974
C   
C   
C     COMPUTE ORTHOGONAL ROTATION MATRIX..  
C     COMPUTE.. MATRIX   (C, S) SO THAT (C, S)(A) = (SQRT(A**2+B**2))   
C                        (-S,C)         (-S,C)(B)   (   0          )    
C     COMPUTE SIG = SQRT(A**2+B**2) 
C        SIG IS COMPUTED LAST TO ALLOW FOR THE POSSIBILITY THAT 
C        SIG MAY BE IN THE SAME LOCATION AS A OR B .
C   
      implicit none
      real a,b,cos,sin,sig,zero,one,xr,yr
      ZERO=0.   
      ONE=1.
      IF (ABS(A).LE.ABS(B)) GO TO 10
      XR=B/A
      YR=SQRT(ONE+XR**2)
      COS=SIGN(ONE/YR,A)
      SIN=COS*XR
      SIG=ABS(A)*YR     
      RETURN
   10 IF (B) 20,30,20   
   20 XR=A/B
      YR=SQRT(ONE+XR**2)
      SIN=SIGN(ONE/YR,B)
      COS=SIN*XR
      SIG=ABS(B)*YR     
      RETURN
   30 SIG=ZERO  
      COS=ZERO  
      SIN=ONE   
      RETURN
      END   


      SUBROUTINE G2    (COS,SIN,X,Y)
C     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1972 DEC 15 
C     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974
C          APPLY THE ROTATION COMPUTED BY G1 TO (X,Y).  
      implicit none
      real cos,sin,x,y,xr
      XR=COS*X+SIN*Y    
      Y=-SIN*X+COS*Y    
      X=XR  
      RETURN
      END   

      function diff(x,y)
      implicit none
      real x,y,diff
      diff=x-y
      return
      end
