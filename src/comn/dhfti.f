      subroutine DHFTI (A,LDA,M1,N1,B,LDB,KB,TAU,KRANK,RNORM,WORK,IP)
c>> 1996-03-30 DHFTI Krogh  Added external statement.
C>> 1994-10-20 DHFTI Krogh  Changes to use M77CON
C>> 1994-04-20 DHFTI CLL Edited to make DP & SP files similar.
c>> 1993-02-09 CLL.  Fixed index in 1st reference to [D/S]NRM2.
C>> 1992-03-13 DHFTI  FTK  Removed implicit statements.
C>> 1987-11-24 DHFTI  Lawson  Initial code.
c--D replaces "?": ?HFTI, ?HTCC, ?HTGEN, ?DOT, ?NRM2
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
c          D1MACH, DHTCC, DHTGEN, DDOT, DNRM2
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
      EXTERNAL D1MACH, DDOT, DNRM2
      INTEGER LDA,M1,N1,LDB,KB,KRANK,IP(N1)
      INTEGER I,II,J,JB,K,KP1,L,LDIAG,LMAX,M,N,NTERMS
      DOUBLE PRECISION D1MACH,DDOT,DNRM2
      DOUBLE PRECISION A(LDA,N1),B(LDB,*),FACTOR,HFAC,ONE
      DOUBLE PRECISION RNORM(KB),SM1,SMALL,TAU,TMP,UPARAM,WORK(N1),ZERO
      logical COMSQR, COL, ROW
      parameter(ONE = 1.0D0, ZERO=0.0D0, FACTOR = 1000.0D0)
      parameter(COL = .true., ROW = .false.)
c     -----------------------------------------------------------------
C
      M = M1
      N = N1
      if( M .lt. 1 .or. N .lt. 0 .or. KB .lt. 0 .or. LDA .lt. M ) then
         call ERMSG('DHFTI',1,0,
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
            call IERM1('DHFTI',2,0,
     *         'Require LDB .ge. 1 when KB .eq. 0', 'LDB', LDB, '.')
            KRANK = 0
            return
         endif
      elseif(LDB .lt. max(M,N)) then
         call IERM1('DHFTI',3,0,
     *      'Require LDB .ge. max(M,N) when KB .ge. 1', 'KB',  KB, ',')
         call IERV1('LDB',LDB, '.')
         KRANK = 0
         return
      endif
c
      if (N .eq. 0) then
         do 10 J = 1, KB
            RNORM(J) = DNRM2(M, B(1,J), 1)
  10     continue
         KRANK = 0
         return
      endif
c                                 Here we have M > 0 and N > 0.
      SMALL = FACTOR * D1MACH(4)
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
         CALL DHTCC (1,N,N+1,M,A(1,N),UPARAM,B,LDB,KB)
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
               WORK(L) = DNRM2(NTERMS, A(J,L), 1)
               IF (WORK(L) .GT. WORK(LMAX)) LMAX = L
   40       continue
            if(WORK(LMAX) .eq. ZERO) then
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
         CALL DHTCC (1,J,J+1,M,A(1,J),UPARAM,A(1,J+1),LDA,N-J)
         CALL DHTCC (2,J,J+1,M,A(1,J),UPARAM,B,LDB,KB)
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
            CALL DHTGEN(1,I,KP1,N,A(I,1),LDA,ROW,WORK(I),A,LDA,I-1,ROW)
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
               SM1 = DDOT(K-I,A(I,I+1),LDA,B(I+1,JB),1)
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
  230          CALL DHTGEN(2,I,KP1,N,A(I,1),LDA,ROW,WORK(I),
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
      
      subroutine AMACH(MODE, I, I1, R1, D1)
c>> 1997-04-16 AMACH  Krogh   Remove blank lines -- was confusing m77con
c>> 1996-03-30 AMACH  Krogh   Added external statement.
c>> 1994-10-26 AMACH  Krogh   Changes to use M77CON
c>> 1994-09-23 AMACH  Snyder  Add VAX G parameters
c>> 1994-06-21 AMACH  Snyder  Compute only round-off and u-flow at first
c>> 1994-05-25 AMACH  Snyder  Added an option to compute at run time.
c>> 1992-04-07 AMACH  Oken    Removed ^Z at EOF (error found by VAX comp
c>> 1992-02-20 AMACH  Snyder  Added Cray-J90 stuff, q.v.
c>> 1990-06-11 AMACH  Snyder  Added Apollo DN-10000 stuff, q.v.
c>> 1990-12-14 AMACH  Lawson  Changed to eliminate ENTRY statements.
c>> 1990-08-21 AMACH  Krogh   No test was getting done for bad machine.
c>> 1990-02-28 AMACH  Krogh   Correct missing DOUBLE PRECISION AMSUB1
c>> 1989-08-14 AMACH  Krogh   Parameterized everything -- Massive change
c>> 1989-03-30 AMACH  Snyder  Correct missing "/" line 921
c>> 1989-01-30 AMACH  Snyder  Incorporate more constants from NETLIB.
C>> 1988-05-19 AMACH  Lawson  Initial code.
c File AMACH.FOR contains user-callable functions I1MACH, D1MACH, and
c DR1MACH , plus second-level subroutines AMACH, AMTEST, and AMSUB1.
c Appropriate lines must be switched between comment and non-comment
c status when this code is moved to a different computer system.
c     These changes can be done with any text editor, however the "c++"
c lines permit automation of the change using the M77CON processor.
c Note that when the M77CON processor activates a line it shifts
c Columns 2-72 to 1-71 and puts a blank in Column 72.  When it inactiv-
c ates a line it shifts Columns 1-71 to 2-72 and puts a C in Column 1.
c     The possible choices using M77CON (don't include parenthetical
c     comments) are:
c      c++ CURRENT HAS SYS = IEEE
c      c++ CURRENT HAS SYS = ALPHA_D3
c      c++ CURRENT HAS SYS = AMDAHL
c      c++ CURRENT HAS SYS = APOLLO_10000
c      c++ CURRENT HAS SYS = BUR1700
c      c++ CURRENT HAS SYS = BUR5700
c      c++ CURRENT HAS SYS = BUR67_7700
c      c++ CURRENT HAS SYS = CDC60_7000
c      c++ CURRENT HAS SYS = CONVEXC_1
c      c++ CURRENT HAS SYS = CRAY1
c      c++ CURRENT HAS SYS = CRAY1_SD (Sngl prec.arith. used for dble.)
c      c++ CURRENT HAS SYS = CRAY1_64 (64 bit integers)
c      c++ CURRENT HAS SYS = CRAY1_SD_64 (64 bit int, SP used for DP)
c      c++ CURRENT HAS SYS = CRAY_T3D
c      c++ CURRENT HAS SYS = CRAY_J90
c      c++ CURRENT HAS SYS = CRAY_J90_SD (Sngl prec. used for dble.)
c      c++ CURRENT HAS SYS = DG_S2000
c      c++ CURRENT HAS SYS = HARRIS220
c      c++ CURRENT HAS SYS = HON600_6000
c      c++ CURRENT HAS SYS = HON_DPS_8_70
c      c++ CURRENT HAS SYS = HP700Q
c      c++ CURRENT HAS SYS = IBM360_370
c      c++ CURRENT HAS SYS = INTERDATA_8_32
c      c++ CURRENT HAS SYS = PDP10_KA
c      c++ CURRENT HAS SYS = PDP10_KB
c      c++ CURRENT HAS SYS = PDP11
c      c++ CURRENT HAS SYS = PRIME50
c      c++ CURRENT HAS SYS = SEQ_BAL_8000
c      c++ CURRENT HAS SYS = UNIVAC
c      c++ CURRENT HAS SYS = VAX
c      c++ CURRENT HAS SYS = VAX_G
c     The current choice is:
c++ CURRENT HAS SYS = IEEE
c
c     One can also select whether floating point constants are created
c     by the compiler or created at run time.  The choices using M77CON
c     are:
c      c++ CURRENT HAS HOW = COMPILER
c      c++ CURRENT HAS HOW = RUN
c     The current choice is:
c++ CURRENT HAS HOW = COMPILER
c
c     If the constants are created at run time, and they fail the run-
c     time check for reasonableness, they are re-created assuming IEEE.
c     If they still fail, the program stops.
c
C  I/O UNIT NUMBERS:
C
C    IM1 = I1MACH( 1) = THE STANDARD INPUT UNIT.
C    IM2 = I1MACH( 2) = THE STANDARD OUTPUT UNIT.
C    IM3 = I1MACH( 3) = THE STANDARD PUNCH UNIT.
C    IM4 = I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
C
C  WORDS:
C
C    IM5 = I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
C    IM6 = I1MACH( 6) = THE NUMBER OF CHARACTERS/INTEGER STORAGE UNIT.
C
C  INTEGERS:
C
C    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM
C
C               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.
C
C    IM7 = I1MACH( 7) = A, THE BASE.
C    IM8 = I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
C    IM9 = I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
C
C  FLOATING-POINT NUMBERS:
C
C    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
C    BASE-B FORM
C
C               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
C               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.
C
C    IM10 = I1MACH(10) = B, THE BASE.
C
C  SINGLE-PRECISION:
C
C    IM11 = I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
C    IM12 = I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
C    IM13 = I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
C
C  DOUBLE-PRECISION:
C
C    IM14 = I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
C    IM15 = I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
C    IM16 = I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
C
C  CONVERSION FROM FUNCTIONAL TO STRUCTURAL FLOATING POINT CONSTANTS
C
C    IM17 = CONSTANT SUCH THAT IM14 + IM17 = ACTUAL NUMBER OF BASE-B
C           DIGITS IN DOUBLE PRECISION, USED FOR CHECKING THAT CORRECT
C           VERSION OF THIS PROGRAM IS INSTALLED.  (SEE DEFINITION OF
C           DM6, AND THE USE OF DM6 IN CALLING AMTEST.)
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF PARAMETER STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.  ALSO, THE VALUES OF
C  IM1 - IM4 SHOULD BE CHECKED FOR CONSISTENCY
C  WITH THE LOCAL OPERATING SYSTEM.
c     -----------------------------------------------------------------
c     Original design and code due to P. A. Fox, A. D. Hall, and
c     N. L. Schryer, Bell Laboratories.  See ACM TOMS, 4,(1978),177-188.
c     Adapted to Univac 1100 by Kris Stewart, JPL, 7/30/81.
c     Adapted for the JPL MATH77 library by C. L. Lawson and F. T. Krogh
c     Sept, 1987.
c     1989-08-14 AMACH  Krogh   Parameterized everything. Major changes.
C     1990 Dec. CLL reorganized code to avoid using ENTRY statements
c     for functions of different types.  Also added save statements.
c     -----------------------------------------------------------------
c     On the first call to this function, tests are done to verify that
c     IM10 and IM14 are not grossly wrong for the host environment.
c     This gives some protection against using the wrong version of this
c     subprogram.
c     -----------------------------------------------------------------
      integer MODE, I, I1
      real R1
      double precision D1, TEST
c
      integer IMACH(17)
      integer IM1, IM2, IM3, IM4, IM5, IM6, IM7, IM8, IM9, IM10, IM11,
     1            IM12, IM13, IM14, IM15, IM16, IM17
c++ Code for HOW=RUN is INACTIVE
C      integer IEEE
C      integer ID1, ID2, ID3, ID4, ID5, ID6, ID7, ID8, ID10, ID11,
C     1   ID12, ID13, ID14, ID15, ID16, ID17
c++ Code for (HOW=RUN) | SYS=IEEE is ACTIVE
      integer IE1, IE2, IE3, IE4, IE5, IE6, IE7, IE8, IE10, IE11,
     1   IE12, IE13, IE14, IE15, IE16, IE17
c++ end
      real             RMACH(5), RM1, RM2, RM3, RM4, RM5,
     1                 RMA, RMB, RBASE
      double precision DMACH(5), DM1, DM2, DM3, DM4, DM5, DM6,
     1                 DMA, DMB, DBASE
      save TEST, IMACH, RMACH, DMACH
C     -----------------------------------------------------------------
C     Machine constants for IEEE standard binary floating-point
c     processors.  This includes PC's and work-stations using the
c     Intel 8087, 80287, 80387, ... processors or the
c     Motorola 68881, 68882, ... processors.
c     Note:  We are setting the "most negative exponent" (IMACH(12) and
c     IMACH(15)) to be the exponent of the smallest normalized number.
c     An IEEE processor actually handles smaller numbers before
c     underflowing, however these "unnormalized" numbers have
c     diminished precision.
c
c++ Code for (HOW=RUN) | SYS=IEEE is ACTIVE
c     Parameters for IEEE when generating at run time:
      PARAMETER (IE1 =5, IE2 =6, IE3 =7, IE4 =6)
      PARAMETER (IE5 =32, IE6 =4, IE7 =2, IE8 =31)
      PARAMETER (IE10 =2, IE11 =24, IE12 =-125, IE13 =128)
      PARAMETER (IE14 =53, IE15 =-1021, IE16 =1024, IE17=0)
c++ Code for SYS = IEEE is ACTIVE
      PARAMETER (IM1 = IE1, IM2 = IE2, IM3 = IE3, IM4 = IE4)
      PARAMETER (IM5 = IE5, IM6 = IE6, IM7 = IE7, IM8 = IE8)
      PARAMETER (IM10 = IE10, IM11 = IE11, IM12 = IE12, IM13 = IE13)
      PARAMETER (IM14 = IE14, IM15 = IE15, IM16 = IE16, IM17 = IE17)
C     -----------------------------------------------------------------
c++ Code for SYS = ALPHA_D3 is INACTIVE
Cc     MACHINE CONSTANTS for the VAX/VMS F and D-3 format for Alpha
Cc
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
C      PARAMETER (IM10 =2, IM11 =24, IM12 =-127, IM13 =127)
C      PARAMETER (IM14 =53, IM15 =-127, IM16 =127, IM17=0)
c++ end
C     -----------------------------------------------------------------
c++ Code for HOW = RUN is INACTIVE
Cc     MACHINE CONSTANTS for the VAX/VMS F and D-3 format for Alpha
Cc
C      PARAMETER (ID1 =5, ID2 =6, ID3 =7, ID4 =6)
C      PARAMETER (ID5 =32, ID6 =4, ID7 =2, ID8 =31)
C      PARAMETER (ID10 =2, ID11 =24, ID12 =-127, ID13 =127)
C      PARAMETER (ID14 =53, ID15 =-127, ID16 =127, ID17=0)
c++ end
C     -----------------------------------------------------------------
c++ Code for SYS = AMDAHL is INACTIVE
CC     MACHINE CONSTANTS FOR AMDAHL MACHINES.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
C      PARAMETER (IM10 =16, IM11 =6, IM12 =-64, IM13 =63)
C      PARAMETER (IM14 =14, IM15 =-64, IM16 =63, IM17=0)
C      -----------------------------------------------------------------
c++ Code for SYS = APOLLO_10000 is INACTIVE
cc     MACHINE CONSTANTS FOR APOLLO DN_10000 MACHINES.
cc     The only difference from IEEE is IM13.  This difference has
cc     nothing to do with the arithmetic or representation used by the
cc     machine.  It is caused by a bug in the compiler:  The right-hand
cc     side of RM2 (below) is apparently evaluated in double precision.
cc     When the compiler is ready to store the resulting value into its
cc     internal data structures, it compares it to an incorrect value
cc     of the overflow limit.  It appears the incorrect value has the
cc     correct exponent, but the fraction is 1.5 instead of 2-2**(-p),
cc     where p is the precision in bits.  You can get the correct result
cc     by changing IM13 to 128, changing RM2 from a parameter to a
cc     variable, and changing the parameter statement that assigns a
cc     value to RM2 into an ordinary assignment statement.
CC
c      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
c      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
c      PARAMETER (IM10 =2, IM11 =24, IM12 =-125, IM13 =127)
c      PARAMETER (IM14 =53, IM15 =-1021, IM16 =1024, IM17 =0)
CC     -----------------------------------------------------------------
c++ Code for SYS = BUR1700 is INACTIVE
CC     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
CC
C      PARAMETER (IM1 =7, IM2 =2, IM3 =2, IM4 =2)
C      PARAMETER (IM5 =36, IM6 =4, IM7 =2, IM8 =33)
C      PARAMETER (IM10 =2, IM11 =24, IM12 =-256, IM13 =255)
C      PARAMETER (IM14 =60, IM15 =-256, IM16 =255, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = BUR5700 is INACTIVE
CC     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =48, IM6 =6, IM7 =2, IM8 =39)
C      PARAMETER (IM10 =8, IM11 =13, IM12 =-50, IM13 =76)
C      PARAMETER (IM14 =26, IM15 =-50, IM16 =76, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = BUR67_7700 is INACTIVE
CC     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =48, IM6 =6, IM7 =2, IM8 =39)
C      PARAMETER (IM10 =8, IM11 =13, IM12 =-50, IM13 =76)
C      PARAMETER (IM14 =26, IM15 =-32754, IM16 =32780, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = CDC60_7000 is INACTIVE
CC     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =60, IM6 =10, IM7 =2, IM8 =48)
C      PARAMETER (IM10 =2, IM11 =47, IM12 =-929, IM13 =1070)
C      PARAMETER (IM14 =94, IM15 =-929, IM16 =1069, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = CONVEXC_1 is INACTIVE
CC     MACHINE CONSTANTS FOR CONVEX C-1.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
C      PARAMETER (IM10 =2, IM11 =24, IM12 =-128, IM13 =127)
C      PARAMETER (IM14 =53, IM15 =-1024, IM16 =1023, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = CRAY1 is INACTIVE
CC     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
C      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =46)
C      PARAMETER (IM10 =2, IM11 =47, IM12 =-8189, IM13 =8190)
C      PARAMETER (IM14 =94, IM15 =-8099, IM16 =8190, IM17=2)
CC     -----------------------------------------------------------------
c++ Code for SYS = CRAY_T3D is INACTIVE
cc     Machine constants for Cray T3D.  IEEE double for both precisions.
c      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
c      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
c      PARAMETER (IM10 =2, IM11 =53, IM12 =-1021, IM13 =1024)
c      PARAMETER (IM14 =53, IM15 =-1021, IM16 =1024, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = CRAY_J90 is INACTIVE
CC     MACHINE CONSTANTS FOR THE CRAY J90
CC     Cray claims the overflow exponent (IM13 and IM16) is 8189, and
CC     the underflow exponent (IM12 and IM15) is -8189, but these values
CC     don't seem to work in cf77:  the underflow limit underflows, and
CC     the overflow limit overflows when using Cray's values.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
C      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =46)
C      PARAMETER (IM10 =2, IM11 =47, IM12 =-8188, IM13 =8189)
C      PARAMETER (IM14 =94, IM15 =-8188, IM16 =8189, IM17=2)
CC     -----------------------------------------------------------------
c++ Code for SYS = CRAY_J90_SD is INACTIVE
CC     MACHINE CONSTANTS FOR THE CRAY J90
CC     Cray claims the overflow exponent (IM13 and IM16) is 8189, and
CC     the underflow exponent (IM12 and IM15) is -8189, but these
CC     values don't seem to work in cf77:  the underflow limit under-
CC     flows, and the overflow limit overflows when using Cray's values.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
C      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =46)
C      PARAMETER (IM10 =2, IM11 =47, IM12 =-8188, IM13 =8189)
C      PARAMETER (IM14 =47, IM15 =-8188, IM16 =8189, IM17=1)
CC     -----------------------------------------------------------------
c++ Code for SYS = CRAY1_SD is INACTIVE
CC     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3, WHEN DOUBLE
CC     PRECISION IS TO USE SINGLE PRECISION ARITHMETIC.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
C      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =46)
C      PARAMETER (IM10 =2, IM11 =47, IM12 =-8189, IM13 =8190)
C      PARAMETER (IM14 =47, IM15 =-8189, IM16 =8190, IM17=1)
CC     -----------------------------------------------------------------
c++ Code for SYS = CRAY1_64 is INACTIVE
CC     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
C      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =63)
C      PARAMETER (IM10 =2, IM11 =47, IM12 =-8189, IM13 =8190)
C      PARAMETER (IM14 =94, IM15 =-8099, IM16 =8190, IM17=2)
CC     -----------------------------------------------------------------
c++ Code for SYS = CRAY1_SD_64 is INACTIVE
CC     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3, WHEN DOUBLE
CC     PRECISION IS TO USE SINGLE PRECISION ARITHMETIC.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
C      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =63)
C      PARAMETER (IM10 =2, IM11 =47, IM12 =-8189, IM13 =8190)
C      PARAMETER (IM14 =47, IM15 =-8189, IM16 =8190, IM17=1)
CC     -----------------------------------------------------------------
c++ Code for SYS = DG_S2000 is INACTIVE
CC     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
CC
C      PARAMETER (IM1 =11, IM2 =12, IM3 =8, IM4 =10)
C      PARAMETER (IM5 =16, IM6 =2, IM7 =2, IM8 =15)
C      PARAMETER (IM10 =16, IM11 =6, IM12 =-64, IM13 =63)
C      PARAMETER (IM14 =14, IM15 =-64, IM16 =63, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = HARRIS220 is INACTIVE
CC     MACHINE CONSTANTS FOR THE HARRIS 220, SLASH 6, SLASH 7.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =0, IM4 =6)
C      PARAMETER (IM5 =24, IM6 =3, IM7 =2, IM8 =23)
C      PARAMETER (IM10 =2, IM11 =23, IM12 =-127, IM13 =127)
C      PARAMETER (IM14 =38, IM15 =-127, IM16 =127, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = HON600_6000 is INACTIVE
CC     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =43, IM4 =6)
C      PARAMETER (IM5 =36, IM6 =6, IM7 =2, IM8 =35)
C      PARAMETER (IM10 =2, IM11 =27, IM12 =-127, IM13 =127)
C      PARAMETER (IM14 =63, IM15 =-127, IM16 =127, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = HON_DPS_8_70 is INACTIVE
CC     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =43, IM4 =6)
C      PARAMETER (IM5 =36, IM6 =4, IM7 =2, IM8 =35)
C      PARAMETER (IM10 =2, IM11 =27, IM12 =-127, IM13 =127)
C      PARAMETER (IM14 =63, IM15 =-127, IM16 =127, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = HP700Q is INACTIVE
cc     Machine constants for HP-700 using the +autodblpad option,
cc     which automatically increases DOUBLE PRECISION to REAL*16, and
cc     REAL to DOUBLE PRECISION.
c      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
c      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
c      PARAMETER (IM10 =2, IM11 =53, IM12 =-1021, IM13 =1024)
c      PARAMETER (IM14 = 113, IM15 = -16381, IM16 = 16384, IM17 = 0)
CC     -----------------------------------------------------------------
c++ Code for SYS = IBM360_370 is INACTIVE
CC     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
CC     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
Cc
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
C      PARAMETER (IM10 =16, IM11 =6, IM12 =-64, IM13 =63)
C      PARAMETER (IM14 =14, IM15 =-64, IM16 =63, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = INTERDATA_8_32 is INACTIVE
CC     MACHINE CONSTANTS FOR THE INTERDATA 8/32
CC     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =6, IM4 =6)
C      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
C      PARAMETER (IM10 =16, IM11 =6, IM12 =-64, IM13 =62)
C      PARAMETER (IM14 =14, IM15 =-64, IM16 =62, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = PDP10_KA is INACTIVE
CC     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =36, IM6 =5, IM7 =2, IM8 =35)
C      PARAMETER (IM10 =2, IM11 =27, IM12 =-128, IM13 =127)
C      PARAMETER (IM14 =54, IM15 =-101, IM16 =127, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = PDP10_KB is INACTIVE
CC     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =36, IM6 =5, IM7 =2, IM8 =35)
C      PARAMETER (IM10 =2, IM11 =27, IM12 =-128, IM13 =127)
C      PARAMETER (IM14 =62, IM15 =-128, IM16 =127, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = PDP11 is INACTIVE
CC     MACHINE CONSTANTS FOR PDP-11 FORTRAN'S SUPPORTING
CC     16-BIT INTEGER ARITHMETIC.
Cc
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =16, IM6 =2, IM7 =2, IM8 =15)
C      PARAMETER (IM10 =2, IM11 =24, IM12 =-127, IM13 =127)
C      PARAMETER (IM14 =56, IM15 =-127, IM16 =127, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = PRIME50 is INACTIVE
CC     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS
CC     WITH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,
CC     SUPPLIED BY IGOR BRAY.
Cc
C      PARAMETER (IM1 =1, IM2 =1, IM3 =2, IM4 =1)
C      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
C      PARAMETER (IM10 =2, IM11 =23, IM12 =-127, IM13 =127)
C      PARAMETER (IM14 =47, IM15 =-32895, IM16 =32637, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = SEQ_BAL_8000 is INACTIVE
CC     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
CC
C      PARAMETER (IM1 =0, IM2 =0, IM3 =7, IM4 =0)
C      PARAMETER (IM5 =32, IM6 =1, IM7 =2, IM8 =31)
C      PARAMETER (IM10 =2, IM11 =24, IM12 =-125, IM13 =128)
C      PARAMETER (IM14 =53, IM15 =-1021, IM16 =1024, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = UNIVAC is INACTIVE
CC     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
CC
CC     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 1
CC     WHICH IS APPROPRIATE FOR THE UNIVAC-FTN SYSTEM.
CC     IF YOU HAVE THE UNIVAC-FOR SYSTEM, SET IT TO 7.
CC     IM6 = 4 for FTN (4 chars per word), 6 for FOR (6 chars per word).
Cc
C      PARAMETER (IM1 =5, IM2 =6, IM3 =1, IM4 =6)
C      PARAMETER (IM5 =36, IM6 =4, IM7 =2, IM8 =35)
C      PARAMETER (IM10 =2, IM11 =27, IM12 =-128, IM13 =127)
C      PARAMETER (IM14 =60, IM15 =-1024, IM16 =1023, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = VAX is INACTIVE
Cc     MACHINE CONSTANTS for the VAX/VMS F and D formats
Cc     and for PDP-11 FORTRAN SUPPORTING 32-BIT INTEGER ARITHMETIC.
Cc
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
C      PARAMETER (IM10 =2, IM11 =24, IM12 =-127, IM13 =127)
C      PARAMETER (IM14 =56, IM15 =-127, IM16 =127, IM17=0)
c++ end
C     -----------------------------------------------------------------
c++ Code for SYS = VAX_G is INACTIVE
Cc     MACHINE CONSTANTS for the VAX/VMS F and G formats
Cc     and for PDP-11 FORTRAN SUPPORTING 32-BIT INTEGER ARITHMETIC.
Cc
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
C      PARAMETER (IM10 =2, IM11 =24, IM12 =-127, IM13 =127)
C      PARAMETER (IM14 =53, IM15 =-1023, IM16 =1023, IM17=0)
c++ end
C     -----------------------------------------------------------------
C
C
C Real parameters
C
C  RM1 = DR1MACH (1) = B**(EMIN-1), The smallest positive number, i.e.,
c                    the underflow limit.
C  RM2 = DR1MACH (2) = B**EMAX*(1 - B**(-T)), The largest number, i.e.,
c                    the overflow limit.
C  RM3 = DR1MACH (3) = B**(-T), The smallest relative spacing, i.e., the
c                    difference between 1.0 and the next smaller number.
C  RM4 = DR1MACH (4) = B**(1-T), The largest relative spacing, i.e., the
c                     difference between 1.0 and the next larger number.
C  RM5 = DR1MACH (5) = LOG10(B).  When B = 2 this value is
c              Log10(2) = 0.30102_99956_63981_19521_37388_94724
C
C Parameter RMA and RMB are selected so that for values of the base =
C 2, 8, 16, 10, RMA has the values 1, 3, 4, 0, and RMB has the values 0,
C 0, 0, 1.  These values are used in computing RM5.
C $$$$ Note that if other bases are to be supported, the calculation of
C $$$$ RMA and RMB will have to be generalized.
C
c++   Code for HOW = COMPILER is ACTIVE
      PARAMETER (IM9 = 2 * (2**(IM8-1) - 1) + 1)
      PARAMETER (RMA = ((IM10 - 10) * (-3 + ((IM10 - 2) * (-77 +
     1    12 * (IM10 - 8))) / 14)) / 24)
      PARAMETER (RMB = ((IM10 - 2) * (IM10 - 8) * (16 - IM10)) / 96)
      PARAMETER (RBASE = IM10)
C
C     Weird subterfuges below are NECESSARY to compute DM1 and DM2 on
C     some systems.  DON'T SIMPLIFY THEM.  We compute RM1 and RM2 using
C     these subterfuges so it will be clear we're computing the REAL
C     and DOUBLE PRECISION characteristics in the same way.
      PARAMETER (RM1 = (RBASE**(IM12/2)) * (RBASE**(IM12-IM12/2-1)))
      PARAMETER (RM2 = RBASE**(IM13-IM11) * ((RBASE**IM11 - RBASE)
     1               + (RBASE - 1.0E0)))
      PARAMETER (RM3 = RBASE**(-IM11))
      PARAMETER (RM4 = RBASE**(1-IM11))
c     PARAMETER (RM5 = RMA*0.30102 99956 63981 19521 37388 94724E0+RMB)
      PARAMETER (RM5 = RMA*0.301029995663981195213738894724E0+RMB)
C
C Double precision parameters -- (Defined like the real ones.)
C
      PARAMETER (DMA = ((IM10 - 10) * (-3 + ((IM10 - 2) * (-77 +
     1    12 * (IM10 - 8))) / 14)) / 24)
      PARAMETER (DMB = ((IM10 - 2) * (IM10 - 8) * (16 - IM10)) / 96)
      PARAMETER (DBASE = IM10)
C
C     Weird subterfuges below are NECESSARY to compute DM1 and DM2 on
C     some systems.  DON'T SIMPLIFY THEM.
      PARAMETER (DM1 = (DBASE**(IM15/2)) * (DBASE**(IM15-IM15/2-1)))
      PARAMETER (DM2 = DBASE**(IM16-IM14) * ((DBASE**IM14 - DBASE)
     1               + (DBASE - 1.0D0)))
      PARAMETER (DM3 = DBASE**(-IM14))
      PARAMETER (DM4 = DBASE**(1-IM14))
c     PARAMETER (DM5 = DMA *
c    1 0.30102 99956 63981 19521 37388 94724 49302 67681 89881 46211 D0
c    2 + DMB)
c
      PARAMETER (DM5 = DMA*
     1 0.30102999566398119521373889472449302676818988146211D0 + DMB)
C DM6 and TEST are used in checking that the correct constants have
C been selected.
      PARAMETER (DM6 = DBASE**(-IM14-IM17))
c++   END
      data TEST / 0.D0 /
C
c     DATA IMACH / IM1, IM2, IM3, IM4, IM5, IM6, IM7, IM8, IM9, IM10,
c    1   IM11, IM12, IM13, IM14, IM15, IM16 /
c     DATA RMACH / RM1, RM2, RM3, RM4, RM5 /
c     DATA DMACH / DM1, DM2, DM3, DM4, DM5 /
C     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      if (TEST .eq. 0.0D0) then
C         IM9 = 2 * (2**(IM8-1) - 1) + 1
         IMACH(1) = IM1
         IMACH(2) = IM2
         IMACH(3) = IM3
         IMACH(4) = IM4
         IMACH(5) = IM5
         IMACH(6) = IM6
         IMACH(7) = IM7
         IMACH(8) = IM8
         IMACH(10) = IM10
         IMACH(11) = IM11
         IMACH(12) = IM12
         IMACH(13) = IM13
         IMACH(14) = IM14
         IMACH(15) = IM15
         IMACH(16) = IM16
         IMACH(17) = IM17
c++   Code for HOW = RUN is INACTIVE
C         IEEE = 0
C100      continue
C      DBASE = IMACH(10)
CC
CC     Weird subterfuge below is NECESSARY to compute DM1 on
CC     some systems.  DON'T SIMPLIFY IT.
C      DM1=(DBASE**(IMACH(15)/2)) * (DBASE**(IMACH(15)-IMACH(15)/2-1))
CC DM6 and TEST are used in checking that the correct constants have
CC been selected.
C      DM6 = DBASE**(-IMACH(14)-IMACH(17))
c++   end
         CALL AMTEST (TEST, DM6)
         if (dm1 .eq. 0.0d0 .or. test .eq. 0.0d0) then
c++   Code for HOW = RUN is INACTIVE
C           if (IEEE .eq. 0) then
C              IEEE = 1
C              IMACH(1) = IE1
C              IMACH(2) = IE2
C              IMACH(3) = IE3
C              IMACH(4) = IE4
C              IMACH(5) = IE5
C              IMACH(6) = IE6
C              IMACH(7) = IE7
C              IMACH(8) = IE8
C              IMACH(10) = IE10
C              IMACH(11) = IE11
C              IMACH(12) = IE12
C              IMACH(13) = IE13
C              IMACH(14) = IE14
C              IMACH(15) = IE15
C              IMACH(16) = IE16
C              IMACH(17) = IE17
C              go to 100
C           end if
C           if (IEEE .eq. 1) then
C              IEEE = 2
C              IMACH(1) = ID1
C              IMACH(2) = ID2
C              IMACH(3) = ID3
C              IMACH(4) = ID4
C              IMACH(5) = ID5
C              IMACH(6) = ID6
C              IMACH(7) = ID7
C              IMACH(8) = ID8
C              IMACH(10) = ID10
C              IMACH(11) = ID11
C              IMACH(12) = ID12
C              IMACH(13) = ID13
C              IMACH(14) = ID14
C              IMACH(15) = ID15
C              IMACH(16) = ID16
C              IMACH(17) = ID17
C              go to 100
C           end if
c++   END
            print*,'AMACH has bad parameters for current environment.'
            stop
         end if
c++   Code for HOW = RUN is INACTIVE
C         IM9 = 2 * (2**(IMACH(8)-1) - 1) + 1
C         RMA = ((IMACH(10) - 10) * (-3 + ((IMACH(10) - 2) * (-77 +
C     1       12 * (IMACH(10) - 8))) / 14)) / 24
C         RMB = ((IMACH(10)-2) * (IMACH(10)-8) * (16-IMACH(10)))/96
C         RBASE = IMACH(10)
CC
CC        Weird subterfuges below are NECESSARY to compute DM1 and DM2
CC        on some systems.  DON'T SIMPLIFY THEM.  We compute RM1 and
CC        RM2 using these subterfuges so it will be clear we're
CC        computing the REAL and DOUBLE PRECISION characteristics in
Cc        the same way.
C         RM1=(RBASE**(IMACH(12)/2))*(RBASE**(IMACH(12)-IMACH(12)/2-1))
C         RM2 = RBASE**(IMACH(13)-IMACH(11))*((RBASE**IMACH(11) - RBASE)
C     1                  + (RBASE - 1.0E0))
C         RM3 = RBASE**(-IMACH(11))
C         RM4 = RBASE**(1-IMACH(11))
Cc        RM5 = RMA*0.30102 99956 63981 19521 37388 94724E0+RMB
C         RM5 = RMA*0.301029995663981195213738894724E0+RMB
CC
CC Double precision parameters -- (Defined like the real ones.)
CC
C         DMA = ((IMACH(10) - 10) * (-3 + ((IMACH(10) - 2) * (-77 +
C     1       12 * (IMACH(10) - 8))) / 14)) / 24
C         DMB = ((IMACH(10)-2) * (IMACH(10)-8) * (16-IMACH(10)))/96
CC
CC        Weird subterfuge below is NECESSARY to compute DM2 on
CC        some systems.  DON'T SIMPLIFY IT.
C         DM2 = DBASE**(IMACH(16)-IMACH(14))*((DBASE**IMACH(14) - DBASE)
C     1                  + (DBASE - 1.0D0))
C         DM3 = DBASE**(-IMACH(14))
C         DM4 = DBASE**(1-IMACH(14))
Cc        DM5 = DMA*0.30102 99956 63981 19521 37388 94724D0+DMB
C         DM5 = DMA*0.301029995663981195213738894724D0+DMB
c++   END
         IMACH(9) = IM9
         RMACH(1) = RM1
         RMACH(2) = RM2
         RMACH(3) = RM3
         RMACH(4) = RM4
         RMACH(5) = RM5
         DMACH(1) = DM1
         DMACH(2) = DM2
         DMACH(3) = DM3
         DMACH(4) = DM4
         DMACH(5) = DM5
      ENDIF
C
      if (MODE .eq. 0) then
         I1=IMACH(I)
      else if (MODE .eq. 1) then
         R1=RMACH(I)
c                                  Here we assume MODE = 2.
      else
         D1=DMACH(I)
      endif
      return
      end
c     ==================================================================
      integer function I1MACH(I)
      integer I, I1
      real R1
      double precision D1
      IF (I .LT. 1  .OR.  I .GT. 16) THEN
         PRINT*,'I1MACH.. Bad argument: I =',I
         STOP 'I1MACH error'
      END IF
      call AMACH (0, I, I1, R1, D1)
      I1MACH = I1
      return
      end
c     ==================================================================
c
      real function DR1MACH (I)
      integer I, I1
      real R1
      double precision D1
      IF (I .lt. 1  .or.  I .gt. 5) THEN
         print*,'DR1MACH .. Bad argument: I = ',I
         stop 'DR1MACH  error'
      END IF
      call AMACH (1, I, I1, R1, D1)
      DR1MACH  = R1
      RETURN
      end
c     ==================================================================
c
      double precision function D1MACH(I)
      integer I, I1
      real R1
      double precision D1
      IF (I .lt. 1  .or.  I .gt. 5) THEN
         print*,'D1MACH.. Bad argument: I = ',I
         stop 'D1MACH error'
      END IF
      call AMACH (2, I, I1, R1, D1)
      D1MACH = D1
      RETURN
      END
c     ==================================================================
c
      SUBROUTINE AMTEST (TEST, D6)
c Verifies that D6 is an appropriate value for DM6.
c Returns TEST = D6 + D6 - 1, .ne. 0 if D6 is an appropriate value for
c DM6, else returns TEST = 0.  The caller uses TEST = 0 as a signal to
c try again with IEEE settings (unless that's already been done).
      external AMSUB1
      DOUBLE PRECISION AMSUB1, D6, TEST
      TEST = AMSUB1(1.D0 + D6)
C
C The comparison with 1.875E0*D6 in the line below is to guard
C against the possibility that TEST is > 0 as a result of rounding
C up in the addition of D6 to 1.
C
      IF ((TEST .eq. 0.D0) .or. (TEST .gt. 1.875D0*D6)) THEN
         TEST = (D6 + D6) + 1.D0
         IF (AMSUB1(TEST) .ne. 0.D0) RETURN
      END IF
      test = 0.0d0
      END
c     ==================================================================
c
      DOUBLE PRECISION FUNCTION AMSUB1 (TEST1)
      DOUBLE PRECISION TEST1
C     Returns the value of TEST1 - 1.
      AMSUB1 = TEST1 - 1.0D0
      RETURN
      END


      DOUBLE PRECISION FUNCTION DDOT(N,X,INCX,Y,INCY)
C>> 1994-11-11 DDOT  Krogh   Declared all vars.
c>> 1994-10-20 DDOT   Krogh  Changes to use M77CON
c>> 1994-04-19 DDOT   Krogh   Minor -- Made code versions line up.
C>> 1985-08-02 DDOT   Lawson  Initial code.
c--D replaces "?": ?DOT
C
C     RETURNS THE DOT PRODUCT OF X AND Y.
C     DDOT = SUM FOR I = 0 TO N-1 OF  X(LX+I*INCX) * Y(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
      INTEGER N, INCX, INCY, IX, IY, I, M, MP1, NS
      DOUBLE PRECISION X(*),Y(*)
      DDOT = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
         DDOT = DDOT + X(IX)*Y(IY)
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
         DDOT = DDOT + X(I)*Y(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         DDOT = DDOT + X(I)*Y(I) + X(I+1)*Y(I+1) +
     $   X(I + 2)*Y(I + 2) + X(I + 3)*Y(I + 3) + X(I + 4)*Y(I + 4)
   50 CONTINUE
      RETURN
C
C         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DDOT = DDOT + X(I)*Y(I)
   70     CONTINUE
      RETURN
      END


      SUBROUTINE DHTCC (MODE,LPIVOT,L1,M,U,  UPARAM,C,LDC,NVC)
c>> 1996-03-30 DHTCC  Krogh   Added external statement.
C>> 1994-11-11 DHTCC  Krogh   Declared all vars.
C>> 1994-10-20 DHTCC  Krogh  Changes to use M77CON
C>> 1987-08-19 DHTCC  Lawson  Initial code.
c--D replaces "?": ?HTCC, ?NRM2
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
c     Subprograms referenced: DNRM2
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
      external DNRM2
      double precision U(*), UPARAM, C(*), DNRM2
      double precision B, FAC, HOLD, VNORM, ONE, ZERO, SUM, BINV
      integer MODE, LPIVOT, L1, M, LDC, NVC
      integer JCBASE, JCPIV, IUL0, J, I
      parameter (ONE = 1.0D0, ZERO = 0.0D0)
C     ------------------------------------------------------------------
      if (0.ge.LPIVOT .or. LPIVOT.ge.L1 .or. L1.gt.M) return
      if( MODE .eq. 1) then
C                            ****** CONSTRUCT THE TRANSFORMATION. ******
         IUL0 = L1 - 1
         if(IUL0 .eq. LPIVOT) then
            VNORM = DNRM2(M-L1+2, U(IUL0), 1)
         else
            HOLD = U(IUL0)
            U(IUL0) = U(LPIVOT)
            VNORM = DNRM2(M-L1+2, U(IUL0), 1)
            U(IUL0) = HOLD
         endif
c
         if (U(LPIVOT) .gt. ZERO) VNORM = -VNORM
         UPARAM = U(LPIVOT)-VNORM
         U(LPIVOT) = VNORM
      endif
C            ****** Apply the transformation  I + U*(U**T)/B  to C. ****
C
      if (NVC .le. 0) return
      B = UPARAM * U(LPIVOT)
C                                 Here B .le. 0.  If B  =  0., return.
      if (B .eq. ZERO) return
      BINV = ONE / B
      JCBASE = 0
c
      do 120 J = 1,NVC
         JCPIV = JCBASE + LPIVOT
         SUM = C(JCPIV) * UPARAM
         do 90 I = L1, M
            SUM = SUM + C(JCBASE+I)*U(I)
   90    continue
         if (SUM .NE. ZERO) then
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


      SUBROUTINE DHTGEN (MODE,LPIVOT,L1,M,U,LDU,COLU,UPARAM,
     *                   C,LDC,NVC,COLC)
c>> 2000-12-06 DHTGEN Krogh For MODE=1 & IUL0=IUPIV, M-L1+2 => M-L1+1
c>> 1996-03-30 DHTGEN Krogh Added external statement.
C>> 1994-11-11 DHTGEN Krogh Declared all vars.
C>> 1994-10-20 DHTGEN Krogh Changes to use M77CON
C>> 1987-08-19 DHTGEN Lawson Initial code.
c--D replaces "?": ?HTGEN, ?AXPY, ?DOT, ?NRM2
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
c     Subprograms referenced: DAXPY, DDOT, DNRM2
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
      external DDOT, DNRM2
      double precision U(*), UPARAM, C(*), DDOT, DNRM2
      double precision B, FAC, HOLD, VNORM, ONE, ZERO, SUM, BINV
      integer MODE, LPIVOT, L1, M, LDU, LDC, NVC
      integer IUPIV, IUL1, IUINC, IUL0
      integer ICE, ICV, I2, I3, INCR, NTERMS, J
      logical COLU, COLC
      parameter (ONE = 1.0D0, ZERO = 0.0D0)
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
            VNORM = DNRM2(M-L1+2, U(IUL0), IUINC)
         else
            HOLD = U(IUL0)
            U(IUL0) = U(IUPIV)
            VNORM = DNRM2(M-L1+2, U(IUL0), IUINC)
            U(IUL0) = HOLD
         endif
c
         if (U(IUPIV) .gt. ZERO) VNORM = -VNORM
         UPARAM = U(IUPIV)-VNORM
         U(IUPIV) = VNORM
      endif
C            ****** Apply the transformation  I + U*(U**T)/B  to C. ****
C
      if (NVC.LE.0) return
      B = UPARAM * U(IUPIV)
c                                 Here B .le. 0.  If B  =  0., return.
      if (B .eq. ZERO) return
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
         SUM = UPARAM * C(I2) + DDOT(NTERMS, U(IUL1),IUINC, C(I3),ICE)
         if (SUM .ne. ZERO) then
            FAC = SUM*BINV
            C(I2) = C(I2) + FAC*UPARAM
            call DAXPY(NTERMS, FAC, U(IUL1),IUINC, C(I3),ICE)
         endif
  120 continue
      return
      end

      DOUBLE PRECISION FUNCTION DNRM2 ( N, X, INCX)
c>> 1998-05-11 DNRM2  Krogh   Minor changes for conversion to C.
c>> 1996-08-29 DNRM2  Krogh   Coded an entirely different algorithm
c>> ....
C>> 1985-08-02 DNRM2  Lawson  Initial code.
C
C     EUCLIDEAN NORM OF THE N-VECTOR STORED IN X() WITH STORAGE
C     INCREMENT INCX .
C     IF    N .LE. 0 RETURN WITH RESULT = 0.
C     IF N .GE. 1 THEN INCX MUST BE .GE. 1
C
C           C.L.LAWSON, 1978 JAN 08
C
C New algorithm avoids underflow as well as overflow, and avoids divides
c as much as possible.  F. Krogh, August 29, 1996.
c
c ************************* Variable Definitions ***********************
c
c Let b be the base for floating point, b**E be the smallest integer
c power of b that overflows, b**e be the smallest integer power of b
c that does not underflow, and d be the number of base b digits in a
c a floating point number.  The descriptions below use b, e, E, and d.
c
c ABIG = b ** ((E+5) / 4)).    This is used as a factor to get the
c   final answer when combining SUMX from big XA, with SUM.
c ASMALL = b ** ((e+d-E+5) / 4)).    This is used as a factor to get the
c   final answer when combining SUMX from small XA, with SUM.
c DBIG = b ** (2*((E+5)/4)).  This is used as a factor to get the
c   final answer when only SUMX formed from big XA's is needed.
c DSMALL = b ** (2*((e+d-E+5)/4)).  This is used as a factor to get the
c   final answer when only SUMX formed from big XA's is needed.
c FBIG = b ** (-2*((E+5)/4)).  This is used as a multiplier when
c   accumulating SUMX for big XA.
c FSMALL = b ** (-2*((e+d-E+5)/4)).  This is used as a multiplier when
c   accumulating SUMX for small XA.
c I      Temporary index.
c ID     Number of base I1MACH(10) digits in Floating point number.
c IEMN   The minimum floating point exponent.
c IEMX   The maximum floating point exponent.
c INCX   Input, the increment (>0) between elements of X.
c N      Input, the number of elements in X.
c NN     N * INCX = last index processed.
c SUM    Place where sum of XA**2 is accumlated.
c SUMX   Place where scaled sum of XA**2 is accumlated.  In first loop
c   this is for XA**2 that would likely underflow, and in second loop
c   this is for XA**2 that would likely overflow.
c TBIG = b ** ((E - d - 1) / 2).  If XA > TBIG, number is "big".  Note
c   that for such XA's, (FBIG * XA) ** 2 should not underflow and the
c   accumlation overflows only if the final result would.
c TSMALL = b ** ((e+1) / 2).  If XA <= TSMALL, number is "small".  Note
c   that for such XA's, (FSMALL * XA)**2 should not underflow and the
c   accumlation should not overlflow.
c X      Input, we are getting the L2 norm of X.
c XA     Contains base of floating point numbers when getting saved
c   parameters.  Later contains abs(X(I)).
C     ------------------------------------------------------------------
c--D replaces "?": ?NRM2
C     ------------------------------------------------------------------
      integer N, INCX
      double precision X(*)
      external I1MACH
      integer I, ID, IEMN, IEMX, I1MACH, NN
      double precision SUM, SUMX, XA
      double precision ABIG,ASMALL,DBIG,DSMALL,FBIG,FSMALL,TBIG,TSMALL
      save ABIG,ASMALL,DBIG,DSMALL,FBIG,FSMALL,TBIG,TSMALL
c Values of these saved parameters for D.P. IEEE arithmetic are:
c   ABIG= .2315841784746324E+78     DBIG= .5363123171977043E+155
c   FBIG= .1864585182800050E-154    TBIG= .9989595361011182E+146
c ASMALL= .4887898181599363E-149  DSMALL= .2389154863368240E-298
c FSMALL= .4185580496821357E+299  TSMALL= .2983336292480080E-153
c Values of these saved parameters for S.P. IEEE arithmetic are:
c   ABIG= .8589935E+10      DBIG= .7378698E+20
c   FBIG= .1355253E-19      TBIG= .2251800E+16
c ASMALL= .1387779E-16    DSMALL= .1925930E-33
c FSMALL= .5192297E+34    TSMALL= .2168404E-18
c
      data ABIG / 0.D0 /
C     ------------------------------------------------------------------
c
c
      if (ABIG .eq. 0.D0) then
C++ Code for (.N. == 'D') is active
         IEMX = I1MACH(16)
         IEMN = I1MACH(15)
         ID = I1MACH(14)
C++ Code for (.N. == 'S') is inactive
C         IEMX = I1MACH(13)
C         IEMN = I1MACH(12)
C         ID = I1MACH(11)
C++ END
         XA = dble(I1MACH(10))
         ABIG = XA ** ((IEMX+5)/4)
         DBIG = ABIG ** 2
         FBIG = 1.D0 / DBIG
         TBIG = XA ** ((IEMX - ID - 1) / 2)
         ASMALL = XA ** ((IEMN + ID - IEMX + 5) / 4)
         DSMALL = ASMALL ** 2
         FSMALL = 1.D0 / DSMALL
         TSMALL = XA ** ((IEMN + 1) / 2)
      end if
      SUM = 0.D0
      if (N .gt. 0) then
         NN = N * INCX
         SUMX = 0.D0
c                      Loop when no big number yet encountered.
         do 100 I = 1, NN, INCX
            XA = abs(X(I))
            if (XA .lt. TSMALL) then
               SUMX = SUMX + (FSMALL * XA) ** 2
            else
               if (XA .gt. TBIG) go to 200
               SUM = SUM + XA**2
            end if
  100    continue
         if (SUM .ne. 0.D0) then
            if (SUMX .ge. 1.D0) then
               if (SUM .lt. 1.D0) then
                  SUM = ASMALL * sqrt(FSMALL*SUM + DSMALL*SUMX)
                  go to 400
               end if
            end if
            SUM = sqrt(SUM)
         else
            SUM = DSMALL * sqrt(SUMX)
         end if
         go to 400
c
  200    SUMX = 0.D0
c                      Loop when we have at least one big number.
         do 300 I = I, NN, INCX
            XA = abs(X(I))
            if (XA .gt. TSMALL) then
               if (XA .gt. TBIG) then
                  SUMX = SUMX + (FBIG * XA) ** 2
               else
                  SUM = SUM + XA**2
               end if
            end if
  300    continue
         if ((SUMX .le. 1.D10) .and. (SUM .ge. 1.D-10)) then
            SUM = ABIG * sqrt(FBIG*SUM + DBIG*SUMX)
         else
            SUM = DBIG * sqrt(SUMX)
         end if
      end if
  400 continue
      DNRM2 = SUM
      return
      end

      SUBROUTINE DAXPY(N,A,X,INCX,Y,INCY)
C>> 1994-11-11 DAXPY  Krogh  Declared all vars.
C>> 1994-10-20 DAXPY  Krogh  Changes to use M77CON
C>> 1985-08-02 DAXPY  Lawson Initial code.
c--D replaces "?": ?AXPY
C
C     OVERWRITE Y WITH A*X + Y.
C     FOR I = 0 TO N-1, REPLACE  Y(LY+I*INCY) WITH A*X(LX+I*INCX) +
C       Y(LY+I*INCY), WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N,
C       AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
C
      INTEGER N, INCX, INCY, IX, IY, I, M, MP1, NS
      DOUBLE PRECISION X(*),Y(*),A
      IF(N.LE.0.OR.A.EQ.0.D0) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
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

