      SUBROUTINE SNNLS(E,NDM,M,N,X,NP,JPVT,Z,NXTRA)
C*********************************************************************
C     This program does a least squares fit of M pieces of data
C     to N unknowns, with possible constraints to a lower limit
C     where E*X=F and E is a M x N matrix with column dimension of NDM,
C     the N+1 column of E contains F, and X is the solution vector of
C     length N. 
C
C     Z is a scratch vector of length N
C
C     ON INPUT:
C     The first N columns of E contain the matrix E
C            and the N+1 column contains F.
C     JPVT(I) is > 0 if parameter I is constrained > 0,
C                    else no constraint
C     M and N should be less than or equal NDM
C     M can be less than N
C
C     ON RETURN:
C     NP is the number of unconstrained parameters
C     X contains the solution.
C     E is upper triangular with the columns permuted like JPVT
C     if K = abs(JPVT(I)), K = input column index and I = new index
C
C     Additional NXTRA columns of data can be added after column N+1
C     of E. These columns will also be transformed.
C*********************************************************************
      implicit none
      INTEGER NDM,M,N,NP,JPVT(*),NXTRA
      REAL E(NDM,*),X(*),Z(*)
      REAL ZERO,EPS,Q
      double precision DSDOT,S,D
      INTEGER NT,I,K,J,JQ,KQ,L,IBGN,IEND,JP1,NDF,MM,NTOT,LAST,JCMP,NX,
     :        NPBGN,ITR,ITRMAX
      PARAMETER (ZERO=0.,EPS=1.E-30)
      LAST=0
      JCMP=0
      NT=N+1
      NTOT=NT+NXTRA
      NDF=N
      IF(NDF.GT.M) NDF=M
      NP=0
      DO 1 J=1,N
          X(J)=ZERO
          IF(JPVT(J).LE.0.AND.NP.LT.NDF) THEN
              NP=NP+1
              JPVT(J)=JPVT(NP)
              JPVT(NP)=J
          ELSE    
              JPVT(J)=J
          ENDIF
   1      CONTINUE
      IBGN=1
      IEND=NP
      NX=IEND
      NPBGN=NP
      MM=M
      ITRMAX=3*NDF
      ITR=0
 200  ITR=ITR+1     
C  begin solution
C********Form Upper Triangle******************************************
        DO 10 J=IBGN,IEND
          JP1 = J+1
C          .. make column norms
          L=MM+1-J
          S=ZERO
          KQ=J
          DO 15 K=J,NX
              JQ=JPVT(K)
              D=DSDOT(L,E(J,JQ),1,E(J,JQ),1)
              IF(D.GT.S) THEN
                S=D
                KQ=K
              ENDIF  
 15           CONTINUE
C          .. swap with current column
          JQ=JPVT(KQ)
          JPVT(KQ)=JPVT(J)
          JPVT(J)=JQ
          S=SQRT(S)
          IF(S.LT.EPS) THEN
              CALL SCOPY(NTOT,ZERO,0,E(J,1),NDM)
              E(J,JQ)=EPS
          ELSE    
              D=E(J,JQ) 
              IF(D.GT.0.) S=-S
              D=D-S
              E(J,JQ)=D
              D=1./(D*S)
              DO 17 K=JP1,N
                  KQ=JPVT(K)
                  Q=D*DSDOT(L,E(J,JQ),1,E(J,KQ),1)
                  CALL SAXPY(L,Q,E(J,JQ),1,E(J,KQ),1)
 17               CONTINUE
              DO 18 K=NT,NTOT
                  Q=D*DSDOT(L,E(J,JQ),1,E(J,K),1)
                  CALL SAXPY(L,Q,E(J,JQ),1,E(J,K),1)
 18               CONTINUE
              E(J,JQ)=S
          ENDIF
          L=L-1
          CALL SCOPY(L,ZERO,0,E(JP1,JQ),1)    
 10       CONTINUE
c        IF(LAST.LT.0 ) write(*,*)'ibgn,iend,last,np,jq,kq,ITR, d='
c        IF(LAST.LT.0 ) write(*,*)ibgn,iend,last,np,jq,kq,itr,d
        IF(LAST.LT.0) RETURN
C*****Compute Solution Vector
        CALL SCOPY(NP,E(1,NT),1,Z,1)
        JQ=0
        D=1.
        DO 19 K=NP,1,-1
          KQ=JPVT(K)
          Z(K)=Z(K)/E(K,KQ)
          Q=-Z(K)
          CALL SAXPY(K-1,Q,E(1,KQ),1,Z,1)
          IF(K.GT.NPBGN) THEN
              IF(Q.GT.ZERO) THEN
                S=X(KQ)/(X(KQ)+Q)
                IF(D.GT.S) THEN
                  JQ=K
                  D=S
                ENDIF
              ENDIF
          ENDIF
 19       CONTINUE
C .....  JQ>0 if there are negative solutions
        IF(JQ.NE.0) THEN
c      write(*,*)'subtracting...np,jq=',np,jq
          JCMP=0
C        .. remove negative solutions from determinable set
          MM=NP
          IBGN=NP
          IEND=NP
          NX=IEND
          K=NP
 20       IF(K.GE.1) THEN  
              I=K-JQ
              KQ=JPVT(K)
              IF(I.NE.0) THEN
                  X(KQ)=X(KQ)+D*( Z(K)-X(KQ) )
                  IF(K.GT.NPBGN) THEN
                      IF(X(KQ).LT.0.) I=0
                  ENDIF
              ENDIF
              IF(I.EQ.0) THEN
                  IF(LAST.EQ.KQ) JCMP=LAST
                  IF(K.LT.NP) THEN
C                  ..swap negative solution to end of determinable set
                      JPVT(K)=JPVT(NP)
                      JPVT(NP)=KQ
                      IBGN=J
                  ENDIF
                  X(NP)=ZERO
                  NP=NP-1
              ENDIF
              K=K-1
              GO TO 20
          ENDIF
          LAST=0
        ELSE
C       .. all current solutions are determinable
          DO 31 K=1,NP
              KQ=JPVT(K)
              X(KQ)=Z(K)
 31           CONTINUE   
c          IF(NP.GE.NDF) write(*,*)'B: itr=',itr
          IF(NP.GE.NDF) RETURN
C*** find largest positive gradient and label with JQ
          JQ=0
          S=0.
          IBGN=NP+1
          MM=M
          L=MM-NP
          DO 33 K=N,IBGN,-1
              KQ=JPVT(K)
              D=DSDOT(L,E(IBGN,KQ),1,E(IBGN,NT),1)
              IF(KQ.NE.JCMP.AND.S.LT.D) THEN
                  JQ=K
                  S=D
              ENDIF
  33          CONTINUE
c      write(*,*)'adding........np,jq=',np,jq
          IF(JQ.EQ.0.OR.ITR.GE.ITRMAX) THEN
C          .. if no positive gradient finish triangularization and quit
              IEND=NDF
              NX=N
              LAST=-1
          ELSE    
C  add largest gradient to determinable set 
              NP=IBGN
              IEND=IBGN
              NX=IEND
              LAST=JPVT(JQ)
              JPVT(JQ)=JPVT(NP)
              JPVT(NP)=LAST
              JCMP=0
          ENDIF    
        ENDIF
      GO TO 200
      END
c

      FUNCTION DSDOT(N,SX,INCX,SY,INCY)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C
      implicit none
      REAL SX(0:*),SY(0:*)
      INTEGER I,INCX,INCY,N,M
      INTEGER IX,IY
      double precision DSDOT
C     EMA SX,SY
C
      DSDOT = 0
      M=N-1
      IF(M.LT.0)RETURN
      IF(INCX.NE.1)GO TO 20
      IF(INCY.NE.1)GO TO 20
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
      DO 30 I = 0,M
        DSDOT = DSDOT + DBLE(SX(I))*SY(I)
   30 CONTINUE
      RETURN
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
   20 IX = 0
      IF(INCX.LT.0)IX = -M*INCX
      IY = 0
      IF(INCY.LT.0)IY = -M*INCY
      DO 10 I = 1,N
        DSDOT = DSDOT + DBLE(SX(IX))*SY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
      END

      FUNCTION SLSFIN(E,NDM,M,N,NPOS,JPVT)
C      returns the Euclidian Norm of the residual vector
C       and the inverse of E (with columns scrambled)
      implicit none
      INTEGER JPVT(*),NDM,M,N,NPOS,I,J,NP1,NT,K,KP1,KQ,JQ,L
      REAL SLSFIN,E(NDM,*),Q,ZERO
      double precision DSDOT
      PARAMETER (ZERO=0.)
      SLSFIN=ZERO
      IF(N.LE.0) RETURN
      NP1=NPOS+1
      NT=N+1
      L=M-NPOS
      SLSFIN=DSQRT( DSDOT(L,E(NP1,NT),1,E(NP1,NT),1))
      NT=MIN(M,N)
c   invert E
      DO 100 K = 1, NPOS
          KQ=JPVT(K)
          E(K,KQ) = 1./E(K,KQ)
          Q = -E(K,KQ)
          CALL SSCAL(K-1,Q,E(1,KQ),1)
          KP1 = K + 1
          DO 80 J = KP1, NPOS
              JQ=JPVT(J)
              Q = E(K,JQ)
              E(K,JQ) = 0.
c              E(J,KQ) = 0.
              CALL SAXPY(K,Q,E(1,KQ),1,E(1,JQ),1)
   80         CONTINUE
  100     CONTINUE
      DO 101 K = NP1,N
          KQ=JPVT(K)
          CALL SCOPY(NT,ZERO,0,E(1,KQ),1)
 101      CONTINUE
c   unscramble if column space available
      IF(N.GT.M) RETURN
      DO 110 I=1,N
c      find new index equivalent to old I
          DO 111 J=I,N
              K=JPVT(J)
              IF(K.EQ.I) GO TO 112
 111          CONTINUE
          J=I
 112      IF(I.NE.J) THEN
c           swap elements so row I contains the proper value
              JPVT(J)=JPVT(I)
              JPVT(I)=I
              CALL SSWAP(N,E(I,1),NDM,E(J,1),NDM)
          ENDIF
 110      CONTINUE
      RETURN
      END
C
      FUNCTION ISAMAX(N,SX,INCX)
C  FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
C
      implicit none
      REAL SX(0:*),SMAX
      INTEGER ISAMAX,N,INCX,IX,I,M
C     EMA SX
C
      ISAMAX = 0
      IF(N.LT.1) RETURN
      M = N-1
      IF(M.NE.0) THEN 
        SMAX = ABS(SX(0))
        IF(INCX.NE.1) THEN
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
          IX = INCX
          DO 10 I = 1,M
             IF(ABS(SX(IX)).GT.SMAX) THEN
                ISAMAX = I
                SMAX = ABS(SX(IX))
             ENDIF
             IX = IX + INCX
   10        CONTINUE
        ELSE
C
C        CODE FOR INCREMENT EQUAL TO 1
C
          DO 30 I = 1,M
             IF(ABS(SX(I)).GT.SMAX) THEN
                 ISAMAX = I 
                 SMAX = ABS(SX(I))
             ENDIF
   30        CONTINUE
        ENDIF
      ENDIF
      ISAMAX = ISAMAX+1
      RETURN
      END

      FUNCTION SASUM(N,SX,INCX)
C  TAKES THE SUM OF THE ABSOLUTE VALUES.
C
      implicit none
      REAL  SASUM,SX(0:*)
      double precision DASUM
      INTEGER I,INCX,N,IX
C     EMA SX
C
      IF(N.LE.0) THEN
          SASUM = 0
          RETURN
      ENDIF
      DASUM = 0
      IF(INCX.NE.1) THEN
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
        IX=0
        DO 10 I = 1,N
          DASUM = DASUM + ABS(SX(IX))
          IX = IX + INCX
   10     CONTINUE
      ELSE
C
C        CODE FOR INCREMENT EQUAL TO 1
C
        IX= N-1
        DO 30 I = 0,IX
          DASUM = DASUM + ABS(SX(I))
   30     CONTINUE
      ENDIF
      SASUM = DASUM
      RETURN
      END
      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C
      implicit none
      double precision da
      REAL  SX(0:*),SY(0:*),SA
      INTEGER I,INCX,INCY,IX,IY,N
C     EMA SX,SY
C
      IF(N.LE.0) RETURN
      IF(ABS(SA).LT.2.e-38) RETURN
      da=sa
      IF(INCX.NE.1.OR.INCY.NE.1) THEN
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
        IX = 0
        IY = 0
        IF(INCX.LT.0) IX = (-N+1)*INCX
        IF(INCY.LT.0) IY = (-N+1)*INCY
        DO 10 I = 1,N
          SY(IY) = SY(IY) + DA*SX(IX)
          IF(ABS(SY(IY)).LT.2.e-38) SY(IY)=0.
          IX = IX + INCX
          IY = IY + INCY
   10     CONTINUE
      ELSE 
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
        IX = N-1
        DO 30 I = 0,IX
           SY(I) = SY(I) + DA*SX(I)
           IF(ABS(SY(I)).LT.2.e-38) SY(I)=0.
   30      CONTINUE
      ENDIF
      RETURN
      END

      SUBROUTINE  SCOPY(N,SX,INCX,SY,INCY)
C COPIES A VECTOR, X, TO A VECTOR, Y.
C
      implicit none
      REAL  SX(0:*),SY(0:*)
      INTEGER I,INCX,INCY,IX,IY,N
C     EMA SX,SY
C
      IF(N.LE.0)RETURN
      IF(INCX.NE.1.OR.INCY.NE.1) THEN
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
        IX = 0
        IY = 0
        IF(INCX.LT.0)IX = (-N+1)*INCX
        IF(INCY.LT.0)IY = (-N+1)*INCY
        DO 10 I = 1,N
          SY(IY) = SX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10     CONTINUE
      ELSE
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
        IX = N-1
        DO 30 I = 0,IX
          SY(I) = SX(I)
   30     CONTINUE
      ENDIF
      RETURN
      END

      FUNCTION SDOT(N,SX,INCX,SY,INCY)
C  FORMS THE DOT PRODUCT OF TWO VECTORS.
C
      implicit none
      REAL SDOT,SX(0:*),SY(0:*)
      double precision DDOT
      INTEGER I,INCX,INCY,IX,IY,N
C     EMA SX,SY
C
      SDOT = 0
      IF(N.LE.0) RETURN
      DDOT=0
      IF(INCX.NE.1.OR.INCY.EQ.1) THEN
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
        IX = 0
        IY = 0
        IF(INCX.LT.0)IX = (-N+1)*INCX
        IF(INCY.LT.0)IY = (-N+1)*INCY
        DO 10 I = 1,N
          DDOT = DDOT + dble(SX(IX))*SY(IY)
          IX = IX + INCX
          IY = IY + INCY
   10     CONTINUE
      ELSE
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
        IX= N-1
        DO 30 I = 0,IX
          DDOT = DDOT + dble(SX(I))*SY(I)
   30     CONTINUE
      ENDIF
      IF(ABS(DDOT).GT.2.e-38) SDOT=DDOT
      RETURN
      END

      SUBROUTINE  SSCAL(N,SA,SX,INCX)
C  SCALES A VECTOR BY A CONSTANT.
C
      implicit none
      REAL  SA,SX(0:*)
      INTEGER I,INCX,N,IX
C     EMA SX
C
      IF(N.LE.0)  RETURN
      IF(SA.EQ.1) RETURN
      IF(INCX.NE.1) THEN
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
          IX=0
          DO 10 I = 1,N
            SX(IX) = SA*SX(IX)
            IX = IX + INCX
   10       CONTINUE
      ELSE
C
C        CODE FOR INCREMENT EQUAL TO 1
C
        IX = N-1
        DO 30 I = 0,IX
          SX(I) = SA*SX(I)
   30     CONTINUE
      ENDIF
      RETURN
      END

      SUBROUTINE  SSWAP (N,SX,INCX,SY,INCY)
C  INTERCHANGES TWO VECTORS.
C
      implicit none
      REAL  SX(0:*),SY(0:*),STEMP
      INTEGER I,INCX,INCY,IX,IY,N
C     EMA SX,SY
C
      IF(N.LE.0)RETURN
      IF(INCX.NE.1.OR.INCY.NE.1) THEN
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
        IX = 0
        IY = 0
        IF(INCX.LT.0)IX = (-N+1)*INCX
        IF(INCY.LT.0)IY = (-N+1)*INCY
        DO 10 I = 1,N
          STEMP = SX(IX)
          SX(IX) = SY(IY)
          SY(IY) = STEMP
          IX = IX + INCX
          IY = IY + INCY
   10     CONTINUE
      ELSE
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
        IX = N-1
        DO 30 I = 0,IX
          STEMP = SX(I)
          SX(I) = SY(I)
          SY(I) = STEMP
   30     CONTINUE
      ENDIF
      RETURN
      END

      SUBROUTINE  SROT (N,SX,INCX,SY,INCY,C,S)
C  APPLIES A PLANE ROTATION.
C
      implicit none
      double precision dc,ds
      REAL  SX(0:*),SY(0:*),C,S,STEMP
      INTEGER INCX,INCY,N,I
      INTEGER IX,IY
C     EMA SX,SY
C
      IF(N.LE.0)RETURN
      dc=c
      ds=s
      IF(INCX.NE.1.OR.INCY.NE.1) THEN
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
        IX = 0
        IY = 0
        IF(INCX.LT.0)IX = (-N+1)*INCX
        IF(INCY.LT.0)IY = (-N+1)*INCY
        DO 10 I = 1,N
          STEMP  = SX(IX)*DC + SY(IY)*DS
          SY(IY) = SY(IY)*DC - SX(IX)*DS
          SX(IX) = STEMP
          IX = IX + INCX
          IY = IY + INCY
   10     CONTINUE
      ELSE
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
        IX = N-1
        DO 30 I = 0,N-1
          STEMP = SX(I)*DC + SY(I)*DS
          SY(I) = SY(I)*DC - SX(I)*DS
          SX(I) = STEMP
   30     CONTINUE
      ENDIF
      RETURN
      END
