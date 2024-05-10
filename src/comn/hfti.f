      subroutine SHFTI (A,LDA,M1,N1,B,LDB,KB,TAU,KRANK,
     & RNORM,WORK,IP)
C>> 1994-10-20 SHFTI Krogh  Changes to use M77CON
C>> 1994-04-20 SHFTI CLL Edited to make DP & SP files similar.
c>> 1993-02-09 CLL.  Fixed index in 1st reference to [D/S]NRM2.
C>> 1992-03-13 SHFTI  FTK  Removed implicit statements.
C>> 1987-11-24 SHFTI  Lawson  Initial code.
c--S replaces "?": ?HFTI, ?HTCC, ?HTGEN, ?DOT, ?NRM2
c
c     ------------------------------------------------------------------
c          This subr solves the least squares problem
c
c                          A * X  ~=~  B
c
c     where A is a given M x N matrix, B is a given M x KB matrix and
c     X is the N x KB solution matrix to be determined.  This includes
c     the usual special case of KB = 1 where B is an M-vector and the
c     solution, X, is an N-vector.
c
c          This subr permits M > N, M = N, or M < N.  This subr
c     determines the "pseudorank", i.e. the estimated rank, of A based
c     on a user-provided tolerance.  If the pseudorank is less than N,
c     the minimal length solution, i.e. the pseudoinverse solution, to
c     the problem is computed.
c
c          Note that this subr can be used to compute the pseudoinverse
c     of a matrix, A.  Set B to the M x M identity matrix and the
c     solution matrix, X, will be the pseudoinverse of A.
c
c          The algorithm is HFTI from the L & H book.  This method does
c     a Householder QR decomposition from the left.  Then if the
c     pseudorank is less than N it does a second Householder QR
c     decomposition from the right.
c
c          The results returned in A(,), RNORM(), and IP() can be used
c     by subroutine SCOV1 or DCOV1 to compute the covariance matrix of
c     the solution vectors.
c     ------------------------------------------------------------------
c                     SUBROUTINE ARGUMENTS
c
c     A(,)     (In/Out)  On input, contains the M x N matrix, A.  Permit
c              M > N, M = N, or M < N.  On return A(,) will contain an
c              upper triangular matrix of order KRANK that can be used
c              by subr _COV2 to compute a covariance matrix when
c              KRANK = N.
c
c     LDA      (In)  The first dimension of the array A(,).
c              Require LDA .ge. M.
c
c     M        (In)  No. of rows of matrices A and B.  Require M .ge. 1.
c
c     N        (In)  No. of columns of matrix A, and rows of matrix X.
c              Require N .ge. 0.
c
c     B(,)     (In/Out)  If KB > 0, the array B(,) must initially
c              contain the right-side matrix, B, having M rows and KB
c              columns.  On return the array B(,) will contain the
c              N x KB solution matrix X.
c              If KB = 0, this subr will not reference the array B(,).
c
c     LDB      (In)  First dimensioning parameter for the array B(,).
c              If KB > 0, require LDB .ge. Max( M, N).
c              If KB = 0, require LDB .ge. 1.
c
c     KB       (In)  No. of columns of the matrices B and X.
c              Require KB .ge. 0.
c              If KB = 0, this subr will not reference the array B(,).
c
c     TAU      (In)  Absolute tolerance parameter provided by user for
c              pseudorank determination.
c
c     KRANK    (Out)  Set by subr to indicate the pseudorank of A.
c              This means that the first KRANK diagonal elements in the
c              the upper triangular factor matrix derived from A each
c              exceed TAU in magnitude.  Either KRANK = Min( M, N), or
c              the the magnitude of the diagonal element in position
c              KRANK + 1 is less than or equal to TAU.
c
c     RNORM()  (Out)  On return, RNORM(J) will contain the euclidean
c              norm of the residual vector for the problem defined by
c              the Jth column vector of the input matrix, B, for
c              J = 1, ..., KB.
c
c     WORK()  (Scratch)  Array used for work space by this subr.
c             Must be of length at least N.
c
c     IP()    (Work/Out)  Integer array of length at least N in which
c              the subr will store column permutation information.
c     -----------------------------------------------------------------
c     Subprograms referenced directly: ERMSG, ERMOR, IERM1, IERV1
c          R1MACH, SHTCC, SHTGEN, SDOT, SNRM2
c     Other subprograms needed: ERFIN
c     -----------------------------------------------------------------
c          This code was originally developed by Charles L. Lawson and
c     Richard J. Hanson at Jet Propulsion Laboratory in 1973.  The
c     original code was described and listed in the book,
c
c                  Solving Least Squares Problems
c                  C. L. Lawson and R. J. Hanson
c                  Prentice-Hall, 1974
c
c     Feb 1985, Mar 1987, June 1987.  C. L. Lawson & S. Y. Chiu, JPL.
c     Adapted code from the Lawson & Hanson book to Fortran 77 for use
c     in the JPL MATH77 library.
c     Changed code to provide oveflow avoidance.
c     Replaced previous scratch arrays H() and G() by WORK().
c     Prefixing subprogram names with S or D for s.p. or d.p. versions.
c     Using generic names for intrinsic functions.
c     Adding calls to BLAS and MATH77 error processing subrs in some
c     program units.
c     ------------------------------------------------------------------
c     1983 Sept 22. CLL added computation of RNORM() for the
c     exceptional case of N = 0.
c     -----------------------------------------------------------------
      implicit none
      INTEGER LDA,M1,N1,LDB,KB,KRANK,IP(N1)
      INTEGER I,II,J,JB,K,KP1,L,LDIAG,LMAX,M,N,NTERMS
      REAL             R1MACH,SDOT,SNRM2
      REAL             A(LDA,N1),B(LDB,*),FACTOR,HFAC,ONE
      REAL             RNORM(KB),SM1,SMALL,TAU,TMP,UPARAM,WORK(N1),ZERO
      logical COMSQR, COL, ROW
      parameter(ONE = 1.0E0, ZERO=0.0E0, FACTOR = 1000.0E0)
      parameter(COL = .true., ROW = .false.)
c     -----------------------------------------------------------------
      HFAC=0.0  ! Avoid compiler warning (may be used uninitialized)
      M = M1
      N = N1
      if( M .lt. 1 .or. N .lt. 0 .or. KB .lt. 0 .or. LDA .lt. M ) then
         call ERMSG('SHFTI',1,0,
     *   'Bad argument values.  Require M .ge. 1, N .ge. 0,', ',')
         call ERMOR('KB .ge. 0, and LDB .ge. M', ',')
         call IERV1('M',  M,   ',')
         call IERV1('N',  N,   ',')
         call IERV1('KB',  KB,   ',')
         call IERV1('LDA',LDA, '.')
         KRANK = 0
         return
      elseif( KB .eq. 0) then
         if(LDB .le. 0) then
            call IERM1('SHFTI',2,0,
     *         'Require LDB .ge. 1 when KB .eq. 0', 'LDB', LDB, '.')
            KRANK = 0
            return
         endif
      elseif(LDB .lt. max(M,N)) then
         call IERM1('SHFTI',3,0,
     *      'Require LDB .ge. max(M,N) when KB .ge. 1', 'KB',  KB, ',')
         call IERV1('LDB',LDB, '.')
         KRANK = 0
         return
      endif
c
      if (N .eq. 0) then
         do 10 J = 1, KB
            RNORM(J) = SNRM2(M, B(1,J), 1)
  10     continue
         KRANK = 0
         return
      endif
c                                 Here we have M > 0 and N > 0.
      SMALL = FACTOR * R1MACH(4)
      LDIAG = MIN(M,N)
C
      DO 80 J = 1,LDIAG
       if(J .eq. N) then
c      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c           Special for J = N.  This case is much simpler than J < N
c           since there are no more columns of A beyond the jth to be
c           considered for interchange or to be triangularized.
c
         IP(N) = N
         CALL SHTCC (1,N,N+1,M,A(1,N),UPARAM,B,LDB,KB)
c      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       else
c      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                                   Here we have J < N.
         IF (J .EQ. 1) THEN
            COMSQR = .TRUE.
         ELSE
c                  Update scaled squared column lengths and set LMAX.
c
            LMAX = J
            DO 20 L = J,N
               WORK(L) = WORK(L) - (HFAC * A(J-1,L))**2
               IF (WORK(L) .GT. WORK(LMAX)) LMAX = L
   20       continue
            COMSQR =  WORK(LMAX) .LE. SMALL
         ENDIF
C
         IF( COMSQR ) THEN
C
C              Compute scaled squared column lengths and set LMAX.
c              Scaling using HFAC protects against overflow of squared
c              numbers.
c
            NTERMS = M - J + 1
            LMAX = J
            DO 40 L = J,N
               WORK(L) = SNRM2(NTERMS, A(J,L), 1)
               IF (WORK(L) .GT. WORK(LMAX)) LMAX = L
   40       continue
            if(abs(WORK(LMAX)) .lt. 1.e-37 ) then
               HFAC = ONE
            else
               HFAC = ONE/WORK(LMAX)
            endif
            do 45 L = J,N
               WORK(L) = (HFAC * WORK(L))**2
   45       continue
         ENDIF
C
C                               DO COLUMN INTERCHANGES IF NEEDED.
C
         IP(J) = LMAX
         IF (IP(J) .NE. J) THEN
            DO 60 I = 1,M
               TMP = A(I,J)
               A(I,J) = A(I,LMAX)
   60          A(I,LMAX) = TMP
            WORK(LMAX) = WORK(J)
         ENDIF
C
C          Compute the J-th transformation and apply it to A and B.
C          Since we treated J = N as a special case we here have J < N
c          so the reference to A(1,J+1) is valid.
c
         CALL SHTCC (1,J,J+1,M,A(1,J),UPARAM,A(1,J+1),LDA,N-J)
         CALL SHTCC (2,J,J+1,M,A(1,J),UPARAM,B,LDB,KB)
c      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       endif
   80 continue
C
C              DETERMINE THE PSEUDORANK, K, USING THE TOLERANCE, TAU.
C
      K = LDIAG
      DO 90 J = 1,LDIAG
         IF (ABS(A(J,J)).LE.TAU) THEN
            K = J - 1
            GO TO 100
         ENDIF
   90 continue
  100 continue
      KP1 = K + 1
C
C                         COMPUTE THE NORMS OF THE RESIDUAL VECTORS.
C
      DO 130 JB = 1,KB
         TMP = ZERO
         DO 120 I = KP1,M
  120       TMP = TMP + B(I,JB)**2
         RNORM(JB) = SQRT(TMP)
  130 continue
C                          Special termination when Pseudorank = 0
      IF (K .EQ. 0) THEN
         DO 155 JB = 1,KB
            DO 150 I = 1,N
              B(I,JB) = ZERO
  150       continue
  155    continue
         KRANK = 0
         RETURN
      ENDIF
C
C               IF THE PSEUDORANK IS LESS THAN N COMPUTE HOUSEHOLDER
C               DECOMPOSITION OF FIRST K ROWS.
C
      IF (K .NE. N) THEN
         DO 170 II = 1,K
            I = KP1-II
            CALL SHTGEN(1,I,KP1,N,A(I,1),LDA,ROW,WORK(I),A,LDA,I-1,ROW)
  170    continue
      ENDIF
C
      DO 260 JB = 1,KB
C
C                        SOLVE THE K BY K TRIANGULAR SYSTEM.
C
         DO 210 L = 1,K
            I = KP1 - L
            IF (I .LT. K) THEN
               SM1 = SDOT(K-I,A(I,I+1),LDA,B(I+1,JB),1)
            ELSE
               SM1 = ZERO
            END IF
            B(I,JB) = (B(I,JB)-SM1) / A(I,I)
  210    continue
C
C     COMPLETE COMPUTATION OF SOLUTION VECTOR.
C    ..
         IF (K .NE. N) THEN
            DO 220 J = KP1,N
  220          B(J,JB) = ZERO
            DO 230 I = 1,K
  230          CALL SHTGEN(2,I,KP1,N,A(I,1),LDA,ROW,WORK(I),
     *                    B(1,JB),LDB,1,COL)
         ENDIF
C                    RE-ORDER THE SOLUTION VECTOR TO COMPENSATE FOR THE
C                    COLUMN INTERCHANGES.
C
         DO 250 J = LDIAG, 1, -1
            IF (IP(J) .NE. J) THEN
               L = IP(J)
               TMP = B(L,JB)
               B(L,JB) = B(J,JB)
               B(J,JB) = TMP
            ENDIF
  250    continue
  260 continue
      KRANK = K
      return
      end


      subroutine SCOV2(A,MDIM,N,IP,VAR,IERR)
C>> 1994-11-11 SCOV2  Krogh   Declared all vars.
C>> 1994-10-20 SCOV2  Krogh  Changes to use M77CON
C>> 1989-10-20 SCOV2  CLL
C>> 1987-11-24 SCOV2  Lawson  Initial code.
c--S replaces "?": ?COV2, ?DOT, ?SWAP
C     Computes upper triangle of covariance matrix,
c     beginning with the triangular matrix, and permutation record, and
c     data error variance computed by _HFTI.
c     Thus, given a matrix, A, represented by the upper triangular
c     matrix in A(), and the permutation record in IP(), and the data
c     error variance, VAR, we compute the upper triangle of the
c     symmetric matrix, C = VAR * ((A**t)*A)**(-1).
C     Adapted from PROG2 in L & H book.
c     ------------------------------------------------------------------
c     Subprograms referenced directly: SDOT, SSWAP, IERM1
c     Other subprograms needed: ERMSG, IERV1, ERFIN, I1MACH
c     ------------------------------------------------------------------
c                 Subroutine Arguments
c
c     A(,) [inout] On entry, contains the upper triangular matrix, A,
c                  in standard, not packed, storage.  This matrix could
c                  have been produced by _HFTI.  On return, contains the
c                  upper triangle of the symmetric unscaled covariance
c                  matrix.  Elements below the diagonal in A(,) will
c                  not be referenced.
c     MDIM [in]    First dimension of A(,).  Require MDIM .ge. N.
c     N [in]       Order of the matrix in A(,)
c     IP() [in]    Permutation record produced by _HFTI.
c     VAR  [in]    Estimate of variance of data error.
c     IERR [out]   Set to 0 if ok.  Set to J > 0 if the (J,J) element
c                  of the given matrix is zero.  In this latter case
c                  the covariance matrix cannot be computed and the
c                  contents of A(,) on return will be meaningless.
C     ------------------------------------------------------------------
c          This code was originally developed by Charles L. Lawson and
c     Richard J. Hanson at Jet Propulsion Laboratory in 1973.  The
c     original code was described and listed in the book,
c
c                  Solving Least Squares Problems
c                  C. L. Lawson and R. J. Hanson
c                  Prentice-Hall, 1974
c
c     Feb, 1985, C. L. Lawson & S. Y. Chan, JPL.  Adapted code from the
c     Lawson & Hanson book to Fortran 77 for use in the JPL MATH77
c     library.
c     Prefixing subprogram names with S or D for s.p. or d.p. versions.
c     Using generic names for intrinsic functions.
c     Adding calls to BLAS and MATH77 error processing subrs in some
c     program units.
C     1989-10-20 CLL Moved integer declaration earlier to avoid warning
c     msg from Cray compiler.
C     ------------------------------------------------------------------
      integer I, IERR, J, N, MDIM, K, IP1, KP1
      integer IP(N)
      real             A(MDIM,N), ONE, TMP, VAR
      real*8 DSDOT
C
      parameter(ONE = 1.0E0)
C     ------------------------------------------------------------------
C     Replace upper triangular matrix U by its inverse, Inv(U)
c
      do 40 J = 1,N
         if(abs(A(J,J)) .lt. 1.e-37 ) then
            call IERM1('SCOV2',1,0,'Jth diagonal elt is zero',
     *      'J',J,'.')
            IERR = J
            return
         end if
         A(J,J) = ONE / A(J,J)
   40 continue

      do 62 I = 1,N-1
         do 60 J = I+1,N
            A(I,J) = -A(J,J) * sngl(DSDOT(J-I,A(I,I),MDIM,A(I,J),1))
   60    continue
   62 continue
C
C     Replace Inv(U) by upper triangle of Trans(Inv(u)) * Inv(U)
c     multiplied by VAR.
c
      do 92 I = 1,N
         do 90 J = I,N
            A(I,J) = VAR * sngl(DSDOT(N-J+1,A(I,J),MDIM,A(J,J),MDIM))
   90    continue
   92 continue
C                                 Permute rows & columns
      do 150 I = N-1, 1, -1
         if (IP(I) .ne. I) then
            K = IP(I)
            TMP = A(I,I)
            A(I,I) = A(K,K)
            A(K,K) = TMP
            if (I .gt. 1) THEN
              CALL SSWAP(I-1,A(1,I),1,A(1,K),1)
            end if
            IP1 = I + 1
            if (IP1 .lt. K) THEN
              CALL SSWAP(K-I-1,A(I,IP1),MDIM,A(IP1,K),1)
            end if
            KP1 = K + 1
            if (KP1 .le. N) THEN
              CALL SSWAP(N-K,A(I,KP1),MDIM,A(K,KP1),MDIM)
            end if
         end if
  150                     continue
      IERR = 0
      end

      SUBROUTINE SSWAP(N,X,INCX,Y,INCY)
C>> 1994-11-11 SSWAP  Krogh   Declared all vars.
C>> 1994-10-20 SSWAP  Krogh  Changes to use M77CON
C>> 1985-08-02 SSWAP  Lawson  Initial code.
c--S replaces "?": ?SWAP
C
C     INTERCHANGE X and Y.
C     FOR I = 0 TO N-1, INTERCHANGE  X(LX+I*INCX) AND Y(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
      INTEGER N, INCX, INCY, IX, IY, I, M, MP1, NS
      REAL             X(*),Y(*),TEMP1,TEMP2,TEMP3
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) THEN
         IF((INCX-1) <0) THEN
           GOTO 5
         ELSEIF ((INCX-1) ==0) THEN
           GOTO 20
         ELSEIF ((INCX-1) >0) THEN
           GOTO 60 
         END IF
      END IF
    5 CONTINUE
C
C       CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        TEMP1 = X(IX)
        X(IX) = Y(IY)
        Y(IY) = TEMP1
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 3.
C
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        TEMP1 = X(I)
        X(I) = Y(I)
        Y(I) = TEMP1
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        TEMP1 = X(I)
        TEMP2 = X(I+1)
        TEMP3 = X(I+2)
        X(I) = Y(I)
        X(I+1) = Y(I+1)
        X(I+2) = Y(I+2)
        Y(I) = TEMP1
        Y(I+1) = TEMP2
        Y(I+2) = TEMP3
   50 CONTINUE
      RETURN
   60 CONTINUE
C
C     CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
      NS = N*INCX
        DO 70 I=1,NS,INCX
        TEMP1 = X(I)
        X(I) = Y(I)
        Y(I) = TEMP1
   70   CONTINUE
      RETURN
      END


      SUBROUTINE SAXPY(N,A,X,INCX,Y,INCY)
C>> 1994-11-11 SAXPY  Krogh   Declared all vars.
C>> 1994-10-20 SAXPY  Krogh  Changes to use M77CON
C>> 1994-10-10 SAXPY  3Krogh  Changes to use M77CON
C>> 1985-08-02 SAXPY  Lawson  Initial code.
c--S replaces "?": ?AXPY
C
C     OVERWRITE Y WITH A*X + Y.
C     FOR I = 0 TO N-1, REPLACE  Y(LY+I*INCY) WITH A*X(LX+I*INCX) +
C       Y(LY+I*INCY), WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N,
C       AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
C
      INTEGER N, INCX, INCY, IX, IY, I, M, MP1, NS
      REAL             X(*),Y(*),A
      IF(N.LE.0. .OR. abs(A).le.0.E0) RETURN
      IF(INCX.EQ.INCY) THEN
        IF((INCX-1) <0) THEN
          GOTO 5
        ELSEIF ((INCX-1) ==0) THEN
          GOTO 20
        ELSEIF ((INCX-1) >0) THEN
          GOTO 60
        ENDIF
      ENDIF
    5 CONTINUE
C
C        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        Y(IY) = Y(IY) + A*X(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        Y(I) = Y(I) + A*X(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        Y(I) = Y(I) + A*X(I)
        Y(I + 1) = Y(I + 1) + A*X(I + 1)
        Y(I + 2) = Y(I + 2) + A*X(I + 2)
        Y(I + 3) = Y(I + 3) + A*X(I + 3)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          Y(I) = A*X(I) + Y(I)
   70     CONTINUE
      RETURN
      END


      SUBROUTINE ERFIN
C     .  Copyright (C) 1989, California Institute of Technology.
C     .  U. S. Government sponsorship under
C     .  NASA contract NAS7-918 is acknowledged.
C>> 1994-11-11 CLL Typing all variables.
C>> 1985-09-23 ERFIN  Lawson  Initial code.
C
      integer idelta, ialpha
      COMMON/M77ERR/IDELTA,IALPHA
      SAVE /M77ERR/
C
 1003 FORMAT(1X,72('$')/' ')
      PRINT 1003
      IF (IALPHA.GE.2) STOP
      RETURN
      END

      REAL             FUNCTION SDOT(N,X,INCX,Y,INCY)
C>> 1994-11-11 SDOT  Krogh   Declared all vars.
c>> 1994-10-20 SDOT   Krogh  Changes to use M77CON
c>> 1994-04-19 SDOT   Krogh   Minor -- Made code versions line up.
C>> 1985-08-02 SDOT   Lawson  Initial code.
c--S replaces "?": ?DOT
C
C     RETURNS THE DOT PRODUCT OF X AND Y.
C     SDOT = SUM FOR I = 0 TO N-1 OF  X(LX+I*INCX) * Y(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
      INTEGER N, INCX, INCY, IX, IY, I, M, MP1, NS
      REAL             X(*),Y(*)
      SDOT = 0.0E0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) THEN
        IF((INCX-1) <0) THEN
          GOTO 5
        ELSEIF ((INCX-1) ==0) THEN
          GOTO 20
        ELSEIF ((INCX-1) >0) THEN
          GOTO 60
        ENDIF
      ENDIF
    5 CONTINUE
C
C         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
         SDOT = SDOT + X(IX)*Y(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         SDOT = SDOT + X(I)*Y(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         SDOT = SDOT + X(I)*Y(I) + X(I+1)*Y(I+1) +
     $   X(I + 2)*Y(I + 2) + X(I + 3)*Y(I + 3) + X(I + 4)*Y(I + 4)
   50 CONTINUE
      RETURN
C
C         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SDOT = SDOT + X(I)*Y(I)
   70     CONTINUE
      RETURN
      END

      REAL*8           FUNCTION DSDOT(N,X,INCX,Y,INCY)
C>> 1994-11-11 SDOT  Krogh   Declared all vars.
c>> 1994-10-20 SDOT   Krogh  Changes to use M77CON
c>> 1994-04-19 SDOT   Krogh   Minor -- Made code versions line up.
C>> 1985-08-02 SDOT   Lawson  Initial code.
c--S replaces "?": ?DOT
C
C     RETURNS THE DOT PRODUCT OF X AND Y.
C     SDOT = SUM FOR I = 0 TO N-1 OF  X(LX+I*INCX) * Y(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
      INTEGER N, INCX, INCY, IX, IY, I, M, MP1, NS
      REAL             X(*),Y(*)
      DSDOT = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) THEN
        IF((INCX-1) <0) THEN
          GOTO 5
        ELSEIF ((INCX-1) ==0) THEN
          GOTO 20
        ELSEIF ((INCX-1) >0) THEN
          GOTO 60
        ENDIF
      ENDIF
    5 CONTINUE
C
C         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
         DSDOT = DSDOT + X(IX)*dble(Y(IY))
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DSDOT = DSDOT + X(I)*dble(Y(I))
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         DSDOT = DSDOT + X(I)*dble(Y(I)) + X(I+1)*dble(Y(I+1)) +
     $   X(I+2)*dble(Y(I+2)) + X(I+3)*dble(Y(I+3)) + X(I+4)*dble(Y(I+4))
   50 CONTINUE
      RETURN
C
C         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DSDOT = DSDOT + X(I)*dble(Y(I))
   70     CONTINUE
      RETURN
      END

      SUBROUTINE SHTGEN (MODE,LPIVOT,L1,M,U,LDU,COLU,
     *                   UPARAM,C,LDC,NVC,COLC)
C>> 1994-11-11 SHTGEN Krogh   Declared all vars.
C>> 1994-10-20 SHTGEN Krogh  Changes to use M77CON
C>> 1987-08-19 SHTGEN Lawson  Initial code.
c--S replaces "?": ?HTGEN, ?AXPY, ?DOT, ?NRM2
c
C        Construction and/or application of a single Householder
C     transformation..     Q = I + U*(U**T)/b
c     where I is the MxM identity matrix, b is a scalar, and U is an
c     M-dimensional Householder vector.
c        All vectors are M-vectors but only the components in positions
c     LPIVOT, and L1 through M, will be referenced.
c        This version, identified by GEN at the end of its name,
c     has the GENerality to handle the options of the U-vector being
c     stored either as a column or row of a matrix, and the vectors in
c     C() may be either column or row vectors.
c     ------------------------------------------------------------------
C                         Subroutine arguments
c
C     MODE  [in]  = 1 OR 2   When MODE = 1 this subr determines the
c           parameters for a Householder transformation and applies
c           the transformation to NVC vectors.  When MODE = 2 this
c           subr applies a previously determined Householder
c           transformation.
C     LPIVOT  [in]  The index of the pivot element.
C     L1,M  [in]  If L1 .le. M  elements LPIVOT and L1 through M will
c           be referenced.  If L1 .gt. M the subroutine returns
c           immediately.  This may be regarded
c           as performing an identity transformation.
C     U()  [inout]  Contains the "pivot" vector.  Typically U() will be
c           a two-dimensional array in the calling program and the pivot
c           vector may be either a column or row in this array.
c           When MODE = 1 this is the vector from which Householder
c           parameters are to be determined.
c           When MODE = 2 this is the result from previous computation
c           with MODE = 1.
c     LDU  [in]  Leading dimensioning parameter for U() in the calling
c           program where U() is a two-dimensional array.  Gives
c           storage spacing between elements in a row of U() when U() is
c           regarded as a two-dimensional array.
C     COLU  [in]  True means the pivot vector is a column of the 2-dim
c           array U().  Thus the successive elements of the pivot vector
c           are at unit storage spacing.
c           False means the pivot vector is a row of the 2-dim array U()
c           Thus the storage spacing between successive elements is LDU.
c     UPARAM  [inout]  Holds a value that supplements the contents
c           of U() to complete the definition of a
c           Householder transformation.  Computed when MODE = 1 and
c           reused when MODE = 2.
c           UPARAM is the pivot component of the Householder U-vector.
C     C()  [inout]   On entry contains a set of NVC M-vectors to which a
c          Householder transformation is to be applied.
c          On exit contains the set of transformed vectors.
c          Typically in the calling program C() will be a 2-dim array
c          with leading dimensioning parameter LDC.
C          These vectors are the columns of an M x NVC matrix in C(,) if
c          COLC = true, and are rows of an NVC x M matrix in C(,) if
c          COLC = false.
C     LDC  [in]  Leading dimension of C(,).  Require LDC .ge. M if
c           COLC = true.  Require LDC .ge. NVC if COLC = false.
C     NVC  [in]  Number of vectors in C(,) to be transformed.
c           If NVC .le. 0 no reference will be made to the array C(,).
c     COLC  [in]  True means the transformations are to be applied to
c           columns of the array C(,).  False means the transformations
c           are to be applied to rows of the array C(,).
c     ------------------------------------------------------------------
c     Subprograms referenced: SAXPY, SDOT, SNRM2
c     ------------------------------------------------------------------
c          This code was originally developed by Charles L. Lawson and
c     Richard J. Hanson at Jet Propulsion Laboratory in 1973.  The
c     original code was described and listed in the book,
c
c                  Solving Least Squares Problems
c                  C. L. Lawson and R. J. Hanson
c                  Prentice-Hall, 1974
c
c     Feb, 1985, C. L. Lawson & S. Y. Chan, JPL.  Adapted code from the
c     Lawson & Hanson book to Fortran 77 for use in the JPL MATH77
c     library.
c     Prefixing subprogram names with S or D for s.p. or d.p. versions.
c     Using generic names for intrinsic functions.
c     Adding calls to BLAS and MATH77 error processing subrs in some
c     program units.
c     July, 1987. CLL.  Changed user interface so method of specifying
c     column/row storage options is more language-independent.
C     ------------------------------------------------------------------
      real             U(*), UPARAM, C(*), SDOT, SNRM2
      real             B, FAC, HOLD, VNORM, ONE, SUM, BINV
      integer MODE, LPIVOT, L1, M, LDU, LDC, NVC
      integer IUPIV, IUL1, IUINC, IUL0
      integer ICE, ICV, I2, I3, INCR, NTERMS, J
      logical COLU, COLC
      parameter (ONE = 1.0E0)
C     ------------------------------------------------------------------
      if (0.ge.LPIVOT .or. LPIVOT.ge.L1 .or. L1.gt.M) return
      if(COLU) then
         IUPIV = LPIVOT
         IUL1 = L1
         IUINC = 1
      else
         IUPIV = 1 + LDU * (LPIVOT-1)
         IUL1 = 1 + LDU * (L1-1)
         IUINC =  LDU
      endif
c
      if( MODE .eq. 1) then
C                            ****** CONSTRUCT THE TRANSFORMATION. ******
         IUL0 = IUL1 - IUINC
         if(IUL0 .eq. IUPIV) then
            VNORM = SNRM2(M-L1+2, U(IUL0), IUINC)
         else
            HOLD = U(IUL0)
            U(IUL0) = U(IUPIV)
            VNORM = SNRM2(M-L1+2, U(IUL0), IUINC)
            U(IUL0) = HOLD
         endif
c
         if (abs(U(IUPIV)) .lt. 1.e-37 ) VNORM = -VNORM
         UPARAM = U(IUPIV)-VNORM
         U(IUPIV) = VNORM
      endif
C            ****** Apply the transformation  I + U*(U**T)/B  to C. ****
C
      if (NVC.LE.0) return
      B = UPARAM * U(IUPIV)
c                                 Here B .le. 0.  If B  =  0., return.
      if (abs(B) .lt. 1.e-37 ) return
      BINV = ONE / B
c                                  I2 = 1 - ICV + ICE*(LPIVOT-1)
c                                  INCR = ICE * (L1-LPIVOT)
      if(COLC) then
         ICE = 1
         ICV = LDC
         I2 = LPIVOT - LDC
         INCR = L1 - LPIVOT
      else
         ICE = LDC
         ICV = 1
         I2 = ICE*(LPIVOT-1)
         INCR = ICE*(L1-LPIVOT)
      endif
c
      NTERMS = M-L1+1
      do 120 J = 1,NVC
         I2 = I2 + ICV
         I3 = I2 + INCR
         SUM = UPARAM * C(I2) + SDOT(NTERMS, U(IUL1),IUINC, C(I3),ICE)
         if (abs(SUM) .gt. 1.e-37) then
            FAC = SUM*BINV
            C(I2) = C(I2) + FAC*UPARAM
            call SAXPY(NTERMS, FAC, U(IUL1),IUINC, C(I3),ICE)
         endif
  120 continue
      return
      end

      SUBROUTINE SHTCC (MODE,LPIVOT,L1,M,U,UPARAM,C,LDC,
     *                  NVC)
C>> 1994-11-11 SHTCC  Krogh   Declared all vars.
C>> 1994-10-20 SHTCC  Krogh  Changes to use M77CON
C>> 1987-08-19 SHTCC  Lawson  Initial code.
c--S replaces "?": ?HTCC, ?NRM2
c
C        Construction and/or application of a single Householder
C     transformation..     Q = I + U*(U**T)/b
c     where I is the MxM identity matrix, b is a scalar, and U is an
c     M-dimensional Householder vector.
c        All vectors are M-vectors but only the components in positions
c     LPIVOT, and L1 through M, will be referenced.
c        This version, identified by CC at the end of its name, is
c     specialized for the Column-Column case, i.e. U() is a vector or
c     a column of a matrix and C() is regarded a containing a set
c     of column vectors to which transformations will be applied.
c     ------------------------------------------------------------------
C                         Subroutine arguments
c
C     MODE  [in]  = 1 OR 2   When MODE = 1 this subr determines the
c           parameters for a Householder transformation and applies
c           the transformation to NVC vectors.  When MODE = 2 this
c           subr applies a previously determined Householder
c           transformation.
C     LPIVOT  [in]  The index of the pivot element.
C     L1,M  [in]  If L1 .le. M  elements LPIVOT and L1 through M will
c           be referenced.  If L1 .gt. M the subroutine returns
c           immediately.  This may be regarded
c           as performing an identity transformation.
C     U()  [inout]  Contains an M-dimensional vector with storage
c           spacing of 1 between elements.
c           When MODE = 1 this is the vector from which Householder
c           parameters are to be determined.
c           When MODE = 2 this is the result from previous computation
c           with MODE = 1.
c     UPARAM  [inout]  Holds a value that supplements the
c           contents of U() to complete the definition of a
c           Householder transformation.  Computed when MODE = 1 and
c           reused when MODE = 2.
c           UPARAM is the pivot component of the Householder U-vector.
C     C()  [inout]  On entry contains a set of NVC M-vectors to which a
c          Householder transformation is to be applied.
c          On exit contains the set of transformed vectors.
C          These vectors are the columns of an M x NVC matrix in C(,).
C     LDC  [in]  Leading dimension of C(,).  Require LDC .ge. M.
C     NVC  [in]  Number of vectors in C(,) to be transformed.
c           If NVC .le. 0 no reference will be made to the array C(,).
c     ------------------------------------------------------------------
c     Subprograms referenced: SNRM2
c     ------------------------------------------------------------------
c          This code was originally developed by Charles L. Lawson and
c     Richard J. Hanson at Jet Propulsion Laboratory in 1973.  The
c     original code was described and listed in the book,
c
c                  Solving Least Squares Problems
c                  C. L. Lawson and R. J. Hanson
c                  Prentice-Hall, 1974
c
c     Feb, 1985, C. L. Lawson & S. Y. Chan, JPL.  Adapted code from the
c     Lawson & Hanson book to Fortran 77 for use in the JPL MATH77
c     library.
c     Prefixing subprogram names with S or D for s.p. or d.p. versions.
c     Using generic names for intrinsic functions.
c     Adding calls to BLAS and MATH77 error processing subrs in some
c     program units.
c     July, 1987. CLL.  Changed user interface so method of specifying
c     column/row storage options is more language-independent.
C     ------------------------------------------------------------------
      real             U(*), UPARAM, C(*), SNRM2
      real             B, FAC, HOLD, VNORM, ONE, SUM, BINV
      integer MODE, LPIVOT, L1, M, LDC, NVC
      integer  IUL0, J, I
      integer*8 jcbase,jcpiv
      parameter (ONE = 1.0E0)
C     ------------------------------------------------------------------
      if (0.ge.LPIVOT .or. LPIVOT.ge.L1 .or. L1.gt.M) return
      if( MODE .eq. 1) then
C                            ****** CONSTRUCT THE TRANSFORMATION. ******
         IUL0 = L1 - 1
         if(IUL0 .eq. LPIVOT) then
            VNORM = SNRM2(M-L1+2, U(IUL0), 1)
         else
            HOLD = U(IUL0)
            U(IUL0) = U(LPIVOT)
            VNORM = SNRM2(M-L1+2, U(IUL0), 1)
            U(IUL0) = HOLD
         endif
c
         if (abs(U(LPIVOT)) .gt. 1.e-37 ) VNORM = -VNORM
         UPARAM = U(LPIVOT)-VNORM
         U(LPIVOT) = VNORM
      endif
C            ****** Apply the transformation  I + U*(U**T)/B  to C. ****
C
      if (NVC .le. 0) return
      B = UPARAM * U(LPIVOT)
C                                 Here B .le. 0.  If B  =  0., return.
      if (abs(B) .lt. 1.e-37) return
      BINV = ONE / B
      JCBASE = 0
c
      do 120 J = 1,NVC
         JCPIV = JCBASE + LPIVOT
         SUM = C(JCPIV) * UPARAM
c         write(*,*)'L1,M,jcbase = ',L1,M,jcbase
         do 90 I = L1, M
            SUM = SUM + C(JCBASE+I)*U(I)
   90    continue
         if (abs(SUM) .gt. 1.e-37) then
            FAC = SUM * BINV
            C(JCPIV) = C(JCPIV) + FAC*UPARAM
            do 110 I =  L1, M
               C(JCBASE+I) = C(JCBASE+I) + FAC*U(I)
  110       continue
         endif
         JCBASE = JCBASE + LDC
  120 continue
      return
      end

      SUBROUTINE ERMSG(SUBNAM,INDIC,LEVEL,MSG,FLAG)
C     .  Copyright (C) 1989, California Institute of Technology.
C     .  U. S. Government sponsorship under
C     .  NASA contract NAS7-918 is acknowledged.
C>> 1994-11-11 ERMSG  Krogh   Declared all vars.
C>> 1992-10-20 ERMSG  WV Snyder  added ERLSET, ERLGET
C>> 1985-09-25 ERMSG  Lawson  Initial code.
C
C     --------------------------------------------------------------
C
C     Four entries: ERMSG, ERMSET, ERLGET, ERLSET
C     ERMSG initiates an error message. This subr also manages the
C     saved value IDELOC and the saved COMMON block M77ERR to
C     control the level of action. This is intended to be the
C     only subr that assigns a value to IALPHA in COMMON.
C     ERMSET resets IDELOC & IDELTA.  ERLGET returns the last value
C     of LEVEL passed to ERMSG.  ERLSET sets the last value of LEVEL.
C     ERLSET and ERLGET may be used together to determine the level
C     of error that occurs during execution of a routine that uses
C     ERMSG.
C
C     --------------------------------------------------------------
C     SUBROUTINE ARGUMENTS
C     --------------------
C     SUBNAM   A name that identifies the subprogram in which
C              the error occurs.
C
C     INDIC    An integer printed as part of the mininal error
C              message. It together with SUBNAM can be used to
C              uniquely identify an error.
C
C     LEVEL    The user sets LEVEL=2,0,or -2 to specify the
C              nominal action to be taken by ERMSG. The
C              subroutine ERMSG contains an internal variable
C              IDELTA, whose nominal value is zero. The
C              subroutine will compute IALPHA = LEVEL + IDELTA
C              and proceed as follows:
C              If (IALPHA.GE.2)        Print message and STOP.
C              If (IALPHA=-1,0,1)      Print message and return.
C              If (IALPHA.LE.-2)       Just RETURN.
C
C     MSG      Message to be printed as part of the diagnostic.
C
C     FLAG     A single character,which when set to '.' will
C              call the subroutine ERFIN and will just RETURN
C              when set to any other character.
C
C     --------------------------------------------------------------
C
C     C.Lawson & S.Chan, JPL, 1983 Nov
C
C     ------------------------------------------------------------------
      INTEGER OLDLEV
      INTEGER IDELOC, LEVEL, IDELTA, IALPHA, INDIC, IDEL
      COMMON/M77ERR/IDELTA,IALPHA
      CHARACTER*(*) SUBNAM,MSG
      CHARACTER*1 FLAG
      SAVE/M77ERR/,IDELOC,OLDLEV
      DATA IDELOC/0/, OLDLEV /0/
      OLDLEV = LEVEL
      IDELTA = IDELOC
      IALPHA = LEVEL + IDELTA
      IF (IALPHA.GE.-1) THEN
c
c            Setting FILE = 'CON' works for MS/DOS systems.
c
c
        WRITE (*,1001) SUBNAM,INDIC
        WRITE (*,*) MSG
        IF (FLAG.EQ.'.') CALL ERFIN
      ENDIF
      RETURN
C
 1001 FORMAT('0',72('$')/' SUBPROGRAM ',A,' REPORTS ERROR NO. ',I4)
C
C
      ENTRY ERMSET(IDEL)
      IDELTA=IDEL
      IDELOC=IDEL
      RETURN
C
C
      ENTRY ERLSET (LEVEL)
      OLDLEV = LEVEL
      RETURN
C
C     
      ENTRY ERLGET (LEVEL)
      LEVEL = OLDLEV
      RETURN
      END

      SUBROUTINE ERMOR(MSG,FLAG)
C     .  Copyright (C) 1989, California Institute of Technology.
C     .  U. S. Government sponsorship under
C     .  NASA contract NAS7-918 is acknowledged.
C>> 1985-09-20 ERMOR  Lawson  Initial code.
C
C     --------------------------------------------------------------
C     SUBROUTINE ARGUMENTS
C     --------------------
C     MSG      Message to be printed as part of the diagnostic.
C
C     FLAG     A single character,which when set to '.' will
C              call the subroutine ERFIN and will just RETURN
C              when set to any other character.
C
C     --------------------------------------------------------------
C
      COMMON/M77ERR/IDELTA,IALPHA
      INTEGER IDELTA,IALPHA
      SAVE /M77ERR/
      CHARACTER*(*) MSG
      CHARACTER*1 FLAG
C
      IF (IALPHA.GE.-1) THEN
        WRITE (*,*) MSG
        IF (FLAG .EQ. '.') CALL ERFIN
      END IF
C
      RETURN
      END

      subroutine IERM1(SUBNAM,INDIC,LEVEL,MSG,LABEL,
     *                 VALUE,FLAG)
C     .  Copyright (C) 1989, California Institute of Technology.
C     .  U. S. Government sponsorship under
C     .  NASA contract NAS7-918 is acknowledged.
C>> 1990-01-18 CLL Added Integer stmt for VALUE.  Typed all variables.
C>> 1985-08-02 IERM1  Lawson  Initial code.
C
      integer INDIC, LEVEL, VALUE
      character*(*) SUBNAM,MSG,LABEL
      character*1 FLAG
      call ERMSG(SUBNAM,INDIC,LEVEL,MSG,',')
      call IERV1(LABEL,VALUE,FLAG)
C
      return
      end

      SUBROUTINE IERV1(LABEL,VALUE,FLAG)
C     .  Copyright (C) 1989, California Institute of Technology.
C     .  U. S. Government sponsorship under
C     .  NASA contract NAS7-918 is acknowledged.
C>> 1985-09-20 IERV1  Lawson  Initial code.
C
C     ------------------------------------------------------------
C     SUBROUTINE ARGUMENTS
C     --------------------
C     LABEL     An identifing name to be printed with VALUE.
C
C     VALUE     A integer to be printed.
C
C     FLAG      See write up for FLAG in ERMSG.
C
C     ------------------------------------------------------------
C
      COMMON/M77ERR/IDELTA,IALPHA
      INTEGER IDELTA,IALPHA,VALUE
      CHARACTER*(*) LABEL
      CHARACTER*1 FLAG
      SAVE /M77ERR/
C
      IF (IALPHA.GE.-1) THEN
        WRITE (*,1002) LABEL,VALUE
        IF (FLAG .EQ. '.') CALL ERFIN
      ENDIF
      RETURN
C
 1002 FORMAT(3X,A,' = ',I5)
      END

      REAL FUNCTION SNRM2(N,DX,INCX)
c     .. Parameters ..
      REAL             ONE         , ZERO
      PARAMETER      ( ONE = 1.0E+0, ZERO = 0.0E+0 )
c     .. Local Scalars ..
      INTEGER          N,INCX
      REAL             ABSXI, NORM, SCALEFACTOR, SSQ, DX(1+(N-1)*INCX)

      IF( N .lt. 1 .or. INCX .lt. 1 ) THEN
         NORM  = ZERO
      ELSE IF( N.eq.1 )THEN
         NORM  = ABS( DX( 1 ) )
      ELSE
         SCALEFACTOR = ZERO
         SSQ   = ONE
         DO IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( abs(DX(IX)) .gt. ZERO ) THEN
               ABSXI = ABS( DX(IX) )
               IF( SCALEFACTOR .lt. ABSXI ) THEN
                  SSQ = ONE + SSQ*( SCALEFACTOR/ABSXI )**2
                  SCALEFACTOR = ABSXI
               ELSE
                  SSQ   = SSQ  + ( ABSXI/SCALEFACTOR )**2
               END IF
            END IF
         END DO
         NORM  = SCALEFACTOR * SQRT( SSQ )
      END IF
      SNRM2 = NORM
      RETURN
      END


      REAL FUNCTION R1MACH(I)
      INTEGER I
C
C  SINGLE-PRECISION MACHINE CONSTANTS
C  R1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C  R1MACH(2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C  R1MACH(3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C  R1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C  R1MACH(5) = LOG10(B)
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
C     needs to be (2) for AUTODOUBLE, HARRIS SLASH 6, ...
      INTEGER SC
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC
      REAL RMACH(5)
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
      INTEGER J, K, L, T3E(3)
      DATA T3E(1) / 9777664 /
      DATA T3E(2) / 5323660 /
      DATA T3E(3) / 46980 /
C  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES,
C  INCLUDING AUTO-DOUBLE COMPILERS.
C  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
C  ON THE NEXT LINE
      DATA SC/0/
C  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
C  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
C          mail netlib@research.bell-labs.com
C          send old1mach from blas
C  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C      DATA RMACH(1) / O402400000000 /
C      DATA RMACH(2) / O376777777777 /
C      DATA RMACH(3) / O714400000000 /
C      DATA RMACH(4) / O716400000000 /
C      DATA RMACH(5) / O776464202324 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C      DATA SMALL(1) /    8388608 /
C      DATA LARGE(1) / 2147483647 /
C      DATA RIGHT(1) /  880803840 /
C      DATA DIVER(1) /  889192448 /
C      DATA LOG10(1) / 1067065499 /, SC/987/
C      DATA RMACH(1) / O00040000000 /
C      DATA RMACH(2) / O17777777777 /
C      DATA RMACH(3) / O06440000000 /
C      DATA RMACH(4) / O06500000000 /
C      DATA RMACH(5) / O07746420233 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C      DATA RMACH(1) / O000400000000 /
C      DATA RMACH(2) / O377777777777 /
C      DATA RMACH(3) / O146400000000 /
C      DATA RMACH(4) / O147400000000 /
C      DATA RMACH(5) / O177464202324 /, SC/987/
C
      IF (SC .NE. 987) THEN
*        *** CHECK FOR AUTODOUBLE ***
         SMALL(2) = 0
         RMACH(1) = 1E13
         IF (SMALL(2) .NE. 0) THEN
*           *** AUTODOUBLED ***
            IF (      SMALL(1) .EQ. 1117925532
     *          .AND. SMALL(2) .EQ. -448790528) THEN
*              *** IEEE BIG ENDIAN ***
               SMALL(1) = 1048576
               SMALL(2) = 0
               LARGE(1) = 2146435071
               LARGE(2) = -1
               RIGHT(1) = 1017118720
               RIGHT(2) = 0
               DIVER(1) = 1018167296
               DIVER(2) = 0
               LOG10(1) = 1070810131
               LOG10(2) = 1352628735
            ELSE IF ( SMALL(2) .EQ. 1117925532
     *          .AND. SMALL(1) .EQ. -448790528) THEN
*              *** IEEE LITTLE ENDIAN ***
               SMALL(2) = 1048576
               SMALL(1) = 0
               LARGE(2) = 2146435071
               LARGE(1) = -1
               RIGHT(2) = 1017118720
               RIGHT(1) = 0
               DIVER(2) = 1018167296
               DIVER(1) = 0
               LOG10(2) = 1070810131
               LOG10(1) = 1352628735
            ELSE IF ( SMALL(1) .EQ. -2065213935
     *          .AND. SMALL(2) .EQ. 10752) THEN
*              *** VAX WITH D_FLOATING ***
               SMALL(1) = 128
               SMALL(2) = 0
               LARGE(1) = -32769
               LARGE(2) = -1
               RIGHT(1) = 9344
               RIGHT(2) = 0
               DIVER(1) = 9472
               DIVER(2) = 0
               LOG10(1) = 546979738
               LOG10(2) = -805796613
            ELSE IF ( SMALL(1) .EQ. 1267827943
     *          .AND. SMALL(2) .EQ. 704643072) THEN
*              *** IBM MAINFRAME ***
               SMALL(1) = 1048576
               SMALL(2) = 0
               LARGE(1) = 2147483647
               LARGE(2) = -1
               RIGHT(1) = 856686592
               RIGHT(2) = 0
               DIVER(1) = 873463808
               DIVER(2) = 0
               LOG10(1) = 1091781651
               LOG10(2) = 1352628735
            ELSE
               WRITE(*,9010)
               STOP 777
               END IF
         ELSE
            RMACH(1) = 1234567.
            IF (SMALL(1) .EQ. 1234613304) THEN
*              *** IEEE ***
               SMALL(1) = 8388608
               LARGE(1) = 2139095039
               RIGHT(1) = 864026624
               DIVER(1) = 872415232
               LOG10(1) = 1050288283
            ELSE IF (SMALL(1) .EQ. -1271379306) THEN
*              *** VAX ***
               SMALL(1) = 128
               LARGE(1) = -32769
               RIGHT(1) = 13440
               DIVER(1) = 13568
               LOG10(1) = 547045274
            ELSE IF (SMALL(1) .EQ. 1175639687) THEN
*              *** IBM MAINFRAME ***
               SMALL(1) = 1048576
               LARGE(1) = 2147483647
               RIGHT(1) = 990904320
               DIVER(1) = 1007681536
               LOG10(1) = 1091781651
            ELSE IF (SMALL(1) .EQ. 1251390520) THEN
*              *** CONVEX C-1 ***
               SMALL(1) = 8388608
               LARGE(1) = 2147483647
               RIGHT(1) = 880803840
               DIVER(1) = 889192448
               LOG10(1) = 1067065499
            ELSE
               DO 10 L = 1, 3
                  J = SMALL(1) / 10000000
                  K = SMALL(1) - 10000000*J
                  IF (K .NE. T3E(L)) GO TO 20
                  SMALL(1) = J
 10               CONTINUE
*              *** CRAY T3E ***
               CALL I1MCRA(SMALL(1), K, 16, 0, 0)
               CALL I1MCRA(LARGE(1), K, 32751, 16777215, 16777215)
               CALL I1MCRA(RIGHT(1), K, 15520, 0, 0)
               CALL I1MCRA(DIVER(1), K, 15536, 0, 0)
               CALL I1MCRA(LOG10(1), K, 16339, 4461392, 10451455)
               GO TO 30
 20            CALL I1MCRA(J, K, 16405, 9876536, 0)
               IF (SMALL(1) .NE. J) THEN
                  WRITE(*,9020)
                  STOP 777
                  END IF
*              *** CRAY 1, XMP, 2, AND 3 ***
               CALL I1MCRA(SMALL(1), K, 8195, 8388608, 1)
               CALL I1MCRA(LARGE(1), K, 24574, 16777215, 16777214)
               CALL I1MCRA(RIGHT(1), K, 16338, 8388608, 0)
               CALL I1MCRA(DIVER(1), K, 16339, 8388608, 0)
               CALL I1MCRA(LOG10(1), K, 16383, 10100890, 8715216)
               END IF
            END IF
 30      SC = 987
         END IF
*     SANITY CHECK
      IF (RMACH(4) .GE. 1.0) STOP 776
      IF (I .LT. 1 .OR. I .GT. 5) THEN
         WRITE(*,*) 'R1MACH(I): I =',I,' is out of bounds.'
         STOP
         END IF
      R1MACH = RMACH(I)
      RETURN
 9010 FORMAT(/' Adjust autodoubled R1MACH by getting data'/
     *' appropriate for your machine from D1MACH.')
 9020 FORMAT(/' Adjust R1MACH by uncommenting data statements'/
     *' appropriate for your machine.')
* /* C source for R1MACH -- remove the * in column 1 */
*#include <stdio.h>
*#include <float.h>
*#include <math.h>
*float r1mach_(long *i)
*{
*	switch(*i){
*	  case 1: return FLT_MIN;
*	  case 2: return FLT_MAX;
*	  case 3: return FLT_EPSILON/FLT_RADIX;
*	  case 4: return FLT_EPSILON;
*	  case 5: return log10((double)FLT_RADIX);
*	  }
*	fprintf(stderr, "invalid argument: r1mach(%ld)\n", *i);
*	exit(1); return 0; /* else complaint of missing return value */
*}
      END
      SUBROUTINE I1MCRA(A, A1, B, C, D)
**** SPECIAL COMPUTATION FOR CRAY MACHINES ****
      INTEGER A, A1, B, C, D
      A1 = 16777216*B + C
      A = 16777216*A1 + D
      END
