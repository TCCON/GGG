      Subroutine CompAbsCH4(path_to_input_files,sgmin,dsig,nsig,
     & T,Ptot,xCH4,xH2O, AbsV, AbsY, AbsW)
C Subroutine computes the absorption coeffcients of CH4 diluted in air 
C in the 2NU3 region
C
        Implicit Double Complex (Z)
        Implicit Double Precision (A-H,O-Y)
        include "../ggg_int_params.f"
        integer lp
C Max number of Spectral points
        Parameter (nSigmx=500000)
        Parameter (Nhitmx=20000,Nnu3mx=76,Nmatrixmx=300000)
C Hitran data
        Dimension Alart(Nhitmx),Srt(Nhitmx),Shift(Nhitmx),Sigr(Nhitmx)
        Character Syml*5,vup*11,vlow*11,
     &  path_to_input_files*(mpath+80)
C Relaxation matrix (for the nu3 band)
        Dimension WT0(Nmatrixmx),BW(Nmatrixmx)
        Dimension Ju1(Nmatrixmx),Ju2(Nmatrixmx)
        Dimension Jl1(Nmatrixmx),Jl2(Nmatrixmx)
        Character SYMnu31(Nmatrixmx)*5,SYMnu32(Nmatrixmx)*5
C Parameters of lines calculated with LM
        Dimension Sig(Nnu3mx),En(Nnu3mx),Jr(Nnu3mx),S2nu3(Nnu3mx)  
        Dimension Delta(Nnu3mx),Gamat(Nnu3mx),Abt(Nnu3mx),Abr(Nnu3mx)
        Dimension Ro0(Nnu3mx),Ro(Nnu3mx),Dip(Nnu3mx),Deg(Nnu3mx)
        Character Supersy*5,Br*1,Abr*1
        Character SY(Nnu3mx)*1,Supersym(Nnu3mx)*5
C Relaxation matrix for the 2nu3 band
        Dimension Zww0(Nnu3mx,Nnu3mx),Zww(Nnu3mx,Nnu3mx)
        Dimension Zvalp(Nnu3mx),Zvecp(Nnu3mx,Nnu3mx)
        Dimension Zivecp(Nnu3mx,Nnu3mx)
        Dimension Zic(Nnu3mx),Zc(Nnu3mx),zWrk(Nnu3mx)
        Dimension EigVR(Nnu3mx),EigVI(Nnu3mx)
        Dimension SSR(Nnu3mx),SSI(Nnu3mx),AlphR(Nnu3mx),AlphI(Nnu3mx)
C Intermediate results
        Dimension Aloh(nSigmx)
        Dimension AbsTV(nSigmx),AbsTW(nSigmx)
        Dimension AbsV(nSigmx),AbsY(nSigmx),AbsW(nSigmx)
        intrinsic dcmplx,dimag
C Results (Absorption Coefficients)
c      Common/CabsV/AbsV(nSigmx)
c      Common/CabsLM/AbsW(nSigmx)
C
C Constants
      Pi=3.14159265d0
      A=1.4388d0
      aMass=16.D-3
      T0=296.d0
      CC=0.1013/(1.38D-23*T)

      if (Nsig.gt.nSigmx) then
      write(6,*)'Set nSigmx at least equal to',Nsig
      stop
      endif
C Doppler width
      om0=6000.d0
      GamD=DSQRT(2.d0*DLOG(2.d0)*8.314/(3d+8**2.0))*
     & DSQRT(T/aMass)*om0
      Cte=DSQRT(DLOG(2.d0))/GamD
      Cte1=Cte/DSQRT(Pi)
C
C Reading of HITRAN 2008 (Frankenberg et al 2008) data base for the considered region********
C Lines that will be accounted by the LM code are skipped here.
C
        lp=len_trim(path_to_input_files)
        Open(60,File=path_to_input_files(:lp)//'06_hit08.par')
c  2745 lines in this linelist. The ones to be LMd have quantum numbers
        ihit=0
68      Read(60,67,end=69)Iso,Sigmar,Srt0,HWairt0,HWselft0,Enr,
     & Bt,S,vup,vlow,Ju,Jl,Syml
67      Format(2x,i1,f12.6,e10.3,10x,f5.4,f5.3,f10.4,
     & f4.2,f8.6,4x,a11,4x,a11,3x,i2,10x,3x,i2,a5)
        if (Sigmar.lt.sgmin.or.sigmar.gt.sgmin+nsig*dsig) go to 68
        idJ=Ju-Jl
        if (abs(idJ).ne.1) m=50 !bidon
        if (idJ.eq.-1) m=-1*Jl
        if (idJ.eq.1) m=Jl+1
        if (Iso.eq.1.and.vup.eq.'0 0 2 0 1F2'.and.vlow.eq.'0 0 0 0 1A1'
     & .and.Syml.ne.'     '.and.abs(m).lt.10) go to 68
        ihit=ihit+1
C Check that dimensions are fine. Stop if not
      if ( ihit .GT. Nhitmx ) Then
      Write(*,2000)                 
      Stop
      End If
2000  Format(//,1x,'************ PROBLEM !!!! ******************',
     /       /,1x,'Arrays in for Line data storage are too small',
     /       /,1x,'raise the value of Nhitmx')
C
      Alart(ihit)=((T0/T)**Bt)*HWairt0
      Qr=QCH4(Iso,T0)/QCH4(Iso,T)
      Srt(ihit)=CC*Srt0*Qr*dexp(-A*Enr*(1./T-1./T0))
     & /(1.d0-dexp(-A*Sigmar/T0))/Sigmar
      Sigr(ihit)=Sigmar
      Shift(ihit)=S
      Goto 68                                  
69    Close (60)
      Nhit=ihit
c      Write(6,*)'Number of lines calculated without LM ',Nhit
C**********************************************************************
      do ihit=1,Nhit
      Sigr(ihit)=Sigr(ihit)+Shift(ihit)*Ptot
      enddo
      do iexp=1,nSig
        Aloh(iexp)=0.d0
      end do

      do ihit=1,Nhit
      do iexp=1,Nsig
        SigC=sgmin+(iexp-1)*DSig
        Cte2=SigC*(1.d0-DEXP(-1.4388d0*SigC/T))

        XX=(Sigr(ihit)-SigC)*Cte
        YY=Alart(ihit)*Ptot*Cte
        Call CPF(XX,YY,WR,WI)

        aa=WR*Cte1*Srt(ihit)
        Aloh(iexp)=Aloh(iexp)+aa*Cte2
      end do
      end do
C**********************************************************************
C Reading of the relaxation matrix of the nu3 band
C
      Open(91,file=path_to_input_files(:lp)//'RMFnu3ai.dat')
      iMatrix=1
92    Read(91,93,end=94)WT0(iMatrix),BW(iMatrix),
     & Ju1(iMatrix),SYMnu31(iMatrix),Jl1(iMatrix),
     & Ju2(iMatrix),SYMnu32(iMatrix),Jl2(iMatrix)
93    format(E15.7,F16.12,4x,i2,2x,a5,1x,i2,5x,i2,2x,a5,1x,i2)
      if (Jl1(iMatrix).eq.Ju1(iMatrix)) go to 92
      if (Jl2(iMatrix).eq.Ju2(iMatrix)) go to 92
      if (Jl1(iMatrix).gt.10.or.Jl2(iMatrix).gt.10) go to 92
      if (WT0(iMatrix).eq.0.d0) go to 92

C Relaxation matrix element at temperature T
      WT0(iMatrix)=WT0(iMatrix)*(T0/T)**BW(iMatrix)
      imatrix=iMatrix+1
      go to 92
94    close(91)
      NMatrix=iMatrix-1
C**********************************************************************
C Reading data for lines calculated with LM
C
      irai=0
      Open(26,file=path_to_input_files(:lp)//'2nu3_hit08.dat',
     & status='old')
c  Only 76 lines in this linelist -- the ones with |m|<10 for which LM is invoked.
19    Read(26,17,end=18)S,St0,E,Bt,Del,Br,Supersy,J,Gat0,GaPine,DelPine
17    Format(3x,f12.6,e11.4,f10.4,1x,f4.2,1x,f8.6,2x,a1,2x,a5,5x,i2
     &,2x,e9.4,3x,f6.4,f8.4)  
      irai=irai+1
      Sig(irai)=S
      Abt(irai)=Bt
      Abr(irai)=Br
      Jr(irai)=J
      Delta(irai)=Del
      En(irai)=E
      S2nu3(irai)=St0
      Sy(irai)=Supersy(3:3)
      Supersym(irai)=Supersy

      If(Supersy(3:3).eq.'F') then
        Deg(irai)=3
      Endif
      If(Supersy(3:3).eq.'E') then
        Deg(irai)=2
      Endif
      If(Supersy(3:3).eq.'A') then
        Deg(irai)=5
      Endif

      Ro(irai)=Deg(irai)*(2.*J+1.)*dexp(-A*E/T)/QCH4(1,T)
      Ro0(irai)=Deg(irai)*(2.*J+1.)*dexp(-A*E/T0)/QCH4(1,T0)
      Qr=QCH4(1,T0)/QCH4(1,T)
      St=CC*S2nu3(irai)/(1.d0-dexp(-A*S/T0))/S
      Dip(irai)=dsqrt(St/Ro0(irai))
      Gamat(irai)=Gat0*(T0/T)**Abt(irai)
      if (GaPine.ne.0.d0) Gamat(irai)=0.985d0*
     & GaPine*(T0/T)**Abt(irai)/9.869d0
      if (DelPine.ne.0.d0) Delta(irai)=DelPine/9.869d0-5.3d-3
      Goto 19
18    Close (26)
      Nraies=irai
c      Write(6,*)'Number of lines calculated with LM: ',Nraies
C**********************************************************************
C Relaxation matrix*********************************
C
c      Write(6,*)'Construction of the relaxation matrix'
C Diagonal elements
      Do 22 irai=1,Nraies
      Zww0(irai,irai)=dcmplx(0.d0,Gamat(irai)*Ptot)
      Sig(irai)=Sig(irai)+Delta(irai)*Ptot
C Off-diagonal elements
      Do 23 iraip=1,Nraies
      If (iraip.eq.irai) goto 23
         Zww0(irai,iraip)=dcmplx(0.d0,0.d0)
         Do 24 imatrix=1,Nmatrix
             If((Supersym(iraip).eq.SYMnu32(iMatrix)).
     &       and.(Supersym(irai).eq.SYMnu31(iMatrix))) then
                  Zww0(irai,iraip)=dcmplx(0.d0,WT0(iMatrix)*Ptot)
                  Goto 25
              Endif
24         Continue
25         Continue
23         Continue
22    Continue

      Do irai=1,Nraies
         Do iraip=1,Nraies
            ZWW(irai,iraip)=Zww0(irai,iraip)
            if (Abr(irai).ne.Abr(iraip)) ZWW(irai,iraip)=0.d0
            if (Jr(irai).ne.Jr(iraip)) ZWW(irai,iraip)=0.d0
         Enddo
         ZWW(irai,irai)=Zww0(irai,irai)+dcmplx(Sig(irai)-om0)
      Enddo
C**********************************************************************
C Matrix diagonalization
C
c      Write(6,*)'Diagonalization'
      Ntest=Nraies
c     Call DEVCCG (Ntest,Zww,Nnu3mx,Zvalp,Zvecp,Nnu3mx)
      Call EigenC(Ntest,Nnu3mx,Zww,zValp,zVecp,zWrk,iFail)
      Do irai=1,Ntest
         EigVR(irai)=DREAL(zValp(irai))
         EigVI(irai)=DIMAG(zValp(irai))
      Enddo

      Do 30 irai=1,Ntest
         zcc=0.d0
         Do 31 iraip=1,Ntest
            zcc=zcc+Ro(iraip)*Dip(iraip)*Zvecp(iraip,irai)
31       Continue
         Zc(irai)=zcc
30    Continue
C inversion
c      Write(6,*)'Inversion'
!      Call DLINCG (Ntest,Zvecp,Nnu3mx,Zivecp,Nnu3mx)
      Call InvMat(Zvecp,Zivecp,Nnu3mx,Ntest,iFail)
      Do irai=1,Ntest
         zicc=0.d0
         Do 33 iraip=1,Ntest
            zicc=zicc+Dip(iraip)*Zivecp(irai,iraip)
33       Continue
         Zc(irai)=Zc(irai)*zicc
      Enddo

      Do irai=1,Ntest
      SSR(irai)=dreal(Zc(irai))
      SSI(irai)=dimag(Zc(irai))
      AlphR(irai)=EigVR(irai)
      AlphI(irai)=EigVI(irai)
      Enddo
C Absorption coeff calculation
c      Write(6,*)'Absorption coeff calculation with Voigt'
      Do 70 iexp=1,nSig
      AbsTV(iexp)=0.d0
      AbsTW(iexp)=0.d0
      AbsV(iexp)=0.d0
      AbsW(iexp)=0.d0
70    Continue

      do 73 irai=1,Ntest
      Do 72 iexp=1,Nsig
         SigC=sgmin+(iexp-1)*dsig
         Cte2=SigC*(1.d0-DEXP(-1.4388*SigC/T))

         XX=(Sig(irai)-SigC)*Cte
         YY=Gamat(irai)*Ptot*Cte
         Call CPF(XX,YY,WR,WI)

         aa=Cte1*xCH4*Ro(irai)*(Dip(irai)**2)*WR
         AbsTV(iexp)=AbsTV(iexp)+aa*Cte2
72    Continue
73    Continue

c       Write(6,*)'Absorption coeff calculation with LM'
      Do 83 irai=1,Ntest
      Do 82 iexp=1,Nsig
         SigC=sgmin+(iexp-1)*dsig
         Cte2=SigC*(1.d0-DEXP(-1.4388*SigC/T))
         XX=(AlphR(irai)+om0-SigC)*Cte
         YY=AlphI(irai)*Cte
         Call CPF(XX,YY,WR,WI)
         aa=Cte1*xCH4*(SSR(irai)*WR-SSI(irai)*WI)
         AbsTW(iexp)=AbsTW(iexp)+aa*Cte2
82    continue
83    continue

      Do 141 iexp=1,Nsig
c         AbsV(iexp)=(AbsTV(iexp)+Aloh(iexp)*xCH4)
c         AbsW(iexp)=(AbsTW(iexp)+Aloh(iexp)*xCH4)
         AbsV(iexp)=Ptot*AbsTV(iexp)              ! Voigt for selected Q & R-branch lines
         AbsY(iexp)=Ptot*(AbsTW(iexp)-AbsTV(iexp)) ! LM-Voigt for selected Q & R-branch lines
         AbsW(iexp)=Ptot*Aloh(iexp)*xCH4          ! Voigt for other non-LM CH4 lines
141   continue
      
      Return
      End Subroutine CompAbsCH4

C******************************************************

      Double precision function Qch4(iso,T)
c ch4 partition function
      implicit DOUBLE PRECISION (a-h,o-z)
      parameter (niso=3,ncoef=4)
      dimension Qcoef(Niso,Ncoef)
c
c...ch4  --  211
      DATA (Qcoef(1,j),j=1,4)/-.17475d+02, .95375d+00,
     +                .39758d-02,-.81837d-06/
c...ch4  --  311
      DATA (Qcoef(2,j),j=1,4)/-.27757d+02, .17264d+01,
     +                .93304d-02,-.48181d-05/
c...ch4  --  212
      DATA (Qcoef(3,j),j=1,4)/-.89810d+03, .44451d+02,
     +                .17474d+00,-.22469d-04/
      if( (t.lt.70.d0) .or. (t.gt.415.d0) )then
      write(6,*)'temp out of range for ch4 partition function'
      stop
        else
        Qch4 = Qcoef(iso,1)
     +       + Qcoef(iso,2)*T
     +       + Qcoef(iso,3)*T*T
     +       + Qcoef(iso,4)*T*T*T
      endif
      RETURN
      END
C------------------------------------------------------

         Subroutine EigenC(nTrue,nDim,zMat,zVal,zVec,zWrk,iFail)
C*********************************************************************
C "EigenC": EIGEN vectors qnd values of Complex matrix
C ..........................................................
C          .   Subroutine compute the Eigenvalues and the  .
C          .   Eigenvectors of a complex operator using a  .
C          .   Hessenberg decomposition and a pivot method .
C          .................................................
C
C Input/Output Parameters of Routine                       
C ---------------------------------
C          nTrue : True dimension of the matrix (order) (Input)    
C           nDim : Number of rows (and columns of the matrix (Input)
C           zMat : Complex matrix to diagonalize (Input).
C           zVec : Complex matrix of the eigenvectors (in columns)
C                  (output).
C           zVal : Complex array of the eigenvalues (output)
C           zWrk : Complex array use for local storage       
C          iFail : should be zero on output. =1 means that search has
C                  not converged (max number of iterations defined in
C                  poutine VecVal
C
C Accessed Files None
C --------------
C
C Called Routines: 'BalanC', 'Hessen', 'VecVal', 'BackTr'
C ---------------
C
C Called By: 'EqvLines' (eventually)
C ---------
C
C Double Precision Complex Version
C
C J.-M. Hartmann, last change 06 March 1997
C*********************************************************************
C
      Implicit None
         Integer nTrue,nDim,iFail
         Integer I,K,L,iCol,iRow
      Double Precision Norm 
      Double Complex zVal,zVec,zMat
      Double Complex zI,zWrk
      Dimension zMat(nDim,nDim)
      Dimension zVec(nDim,nDim)
      Dimension zVal(nDim),zWrk(nDim)
      Intrinsic dconjg
C
      zI=(0.d0,1.d0)
c- balance matrix     
      Call BalanC(nTrue,nDim,zMat,K,L,zWrk)
      Do 1 I=1,nTrue
      zWrk(I)=zI*zWrk(I)
 1    Continue
c- Put under Hessenberg form        
      Call Hessen(K,L,nTrue,nDim,zMat,zWrk)
c- Search Eigenvectors and Eigenvalues  
      Call VecVal(K,L,nTrue,nDim,zMat,zVec,zVal,zWrk,iFail)    
         If( iFail.NE.0 )Return
c- Backtransformation              
      Call BackTr(nTrue,nDim,zVec,K,L,nTrue,zWrk)
C- Normalize
      Do 2 iCol=1,nTrue
        Norm=0.d0
        Do 20 iRow=1,nTrue
        Norm=Norm+zVec(iRow,iCol)*DCONJG(zVec(iRow,iCol))
  20    Continue
      Norm=DSQRT(Norm)
      Do 21 iRow=1,nTrue
      zVec(iRow,iCol)=zVec(iRow,iCol)/Norm
  21  Continue
 2    Continue
      Return
      End
      Subroutine BackTr(nTrue,nDim,zVec,K,L,M,zWrk)
C*********************************************************************
C "BackTr": BACK Transform                     
C ..........................................................
C          .................................................
C               Part of a Package for Diaginalization of
C                        Complex Operators
C          .................................................
C
C
C Accessed Files None
C --------------
C
C Called Routines: None
C ---------------
C
C Called By: 'EigenC' 
C ---------
C
C Double Precision Complex Version
C
C J.-M. Hartmann, last change 07 March 1997
C*********************************************************************
      Implicit None
        Integer nTrue,nDim
        Integer K,L,M,I,J,II,KK
      Double Precision S
         Double Complex zVec,zWrk,zDum
      Dimension zVec(nDim,nDim),zWrk(nDim)
      Intrinsic dimag
C
      If( L.NE.K )Then   
      Do 1 I=K,L
      S=DIMAG(zWrk(I))
      Do 10 J=1,M
      zVec(I,J)=zVec(I,J)*S
10    Continue
 1    Continue
      End If
C
      Do 2 II=1,nTrue
      I=II
      If( I.GE.K.AND.I.LE.L )GoTo 2
      If( I.LT.K )I=K-II
      KK=DIMAG(zWrk(I))
      If( KK.NE.I )Then   
      Do 20 J=1,M
      zDum=zVec(I,J)
      zVec(I,J)=zVec(KK,J)
      zVec(KK,J)=zDum
20    Continue
      End If
2     Continue
      Return
      End
      Subroutine BalanC(nTrue,nDim,zMat,K,L,zWrk)
C*********************************************************************
C "BalanC": BALANCe the matrix                 
C ..........................................................
C          .................................................
C               Part of a Package for Diaginalization of
C                        Complex Operators
C          .................................................
C
C The parameter MachBa is the MACHine BAse representation. It depends
C on the computer but is TWO on most machines
C
C Accessed Files None
C --------------
C
C Called Routines: None
C ---------------
C
C Called By: 'EigenC' 
C ---------
C
C Double Precision Complex Version
C
C J.-M. Hartmann, last change 07 March 1997
C*********************************************************************
      Implicit None
        Integer nTrue,nDim
        Integer K,L,M,J,I,IEXC,L1,JJ
      Double Precision A,R,Zero,One,Pt95,MachBa,MachM1,B2,B2M1,G,F,S
         Double Complex zMat,zWrk
         Double Complex zZero,zOne,zDum
      Dimension zMat(nDim,nDim),zWrk(nDim)
      Logical NoConv
C
      Zero=0.d0
      One=1.d0
      Pt95=0.95d0
      zZero=(0.D+0,0.D+0)
      zOne=(1.D+0,0.D+0)
C MachBa is the BAse of the MACHibe floating point representation
      MachBa=2.d0
C
      B2=MachBa*MachBa
      MachM1=One/MachBa
      B2M1=MachM1*MachM1
      K=1
      L=nTrue
      GoTo 30
C
  5   Continue
      zWrk(M)=J
      If( J.EQ.M )GoTo 20
      Do 10 I=1,L
      zDum=zMat(I,J)
      zMat(I,J)=zMat(I,M)
      zMat(I,M)=zDum
10    Continue
      Do 15 I=K,nTrue
      zDum=zMat(J,I)
      zMat(J,I)=zMat(M,I)
      zMat(M,I)=zDum
15    Continue
20    Continue
      GoTo (25,45),IEXC
C
25    Continue
      If( L.EQ.1 )Return
      L=L-1
30    Continue
      L1=L+1
      Do 40 JJ=1,L
      J=L1-JJ
      Do 35 I=1,L
      If( I.EQ.J )GoTo 35
      If( CDABS(zMat(I,J)).NE.Zero )GoTo 40
35    Continue
      M=L
      IEXC=1
      GoTo 5
40    Continue
      GoTo 50
C
45    Continue
      K=K+1
50    Continue
      Do 60 J=K,L
      Do 55 I=K,L
      If( I.EQ.J )GoTo 55
      If( CDABS(zMat(I,J)).NE.Zero )GoTo 60
55    Continue
      M=K
      IEXC=2
      GoTo 5
60    Continue
C
      Do 65 I=K,L
      zWrk(I)=One
65    Continue
70    Continue   
      NoConv=.FALSE.
      Do 110 I=K,L
      A=Zero
      R=Zero
      Do 75 J=K,L
      If( J.EQ.I )GoTo 75
      A=A+CDABS(zMat(J,I))
      R=R+CDABS(zMat(I,J))
75    Continue
      G=R*MachM1
      F=One
      S=A+R
80    Continue
      If( A.GE.G )GoTo 85
      F=F*MachBa
      A=A*B2
      GoTo 80
85    Continue
      G=R*MachBa
90    Continue   
      If( A.LT.G )GoTo 95
      F=F*MachM1
      A=A*B2M1
      GoTo 90
C
95    Continue   
      If( (A+R)/F.GE.Pt95*S )GoTo 110
      G=One/F
      zWrk(I)=zWrk(I)*F
      NoConv=.TRUE.
      Do 100 J=K,nTrue
      zMat(I,J)=zMat(I,J)*G
100   Continue   
      Do 105 J=1,L
      zMat(J,I)=zMat(J,I)*F
105   Continue   
110   Continue
      If( NoConv )GoTo 70
      Return
      End
      Subroutine Hessen(K,L,nTrue,nDim,zMat,zWrk)
C*********************************************************************
C "Hessen": make HESSENberg decomposition      
C ..........................................................
C          .................................................
C               Part of a Package for Diaginalization of
C                        Complex Operators
C          .................................................
C
C
C Accessed Files None
C --------------
C
C Called Routines: None
C ---------------
C
C Called By: 'EigenC' 
C ---------
C
C Double Precision Complex Version
C
C J.-M. Hartmann, last change 07 March 1997
C*********************************************************************
C
      Implicit None
        Integer nTrue,nDim
        Integer K,L,LA,KP1,M,I,J,MM1,MP1
      Double Complex zMat,zWrk
      Double Complex zZero,zX,zY,zDum
      Dimension zMat(nDim,nDim),zWrk(nDim)
C
      zZero=(0.d0,0.d0)
      LA=L-1
      KP1=K+1
      If( LA.LT.KP1 )Return
      Do 1 M=KP1,LA
      I=M
      zX=zZero
      Do 10 J=M,L
      If( CDABS(zMat(J,M-1)).LE.CDABS(zX) )GoTo 10
      zX=zMat(J,M-1)
      I=J
 10   Continue
      zWrk(M)=I+zWrk(M)
      If( I.NE.M )Then   
      MM1=M-1
      Do 11 J=MM1,nTrue
      zDum=zMat(I,J)
      zMat(I,J)=zMat(M,J)
      zMat(M,J)=zDum
 11   Continue
      Do 12 J=1,L
      zDum=zMat(J,I)
      zMat(J,I)=zMat(J,M)
      zMat(J,M)=zDum
 12   Continue
      End If
      If( zX.EQ.zZero )GoTo 1
      MP1=M+1
      Do 13 I=MP1,L
      zY=zMat(I,M-1)
      If( zY.NE.zZero )Then   
      zY=zY/zX
      zMat(I,M-1)=zY
      Do 130 J=M,nTrue
      zMat(I,J)=zMat(I,J)-zY*zMat(M,J)
130   Continue
      Do 131 J=1,L
      zMat(J,M)=zMat(J,M)+zY*zMat(J,I)
131   Continue
      End If
13    Continue
 1    Continue
      End
      Subroutine VecVal(K,L,nTrue,nDim,zMat,zVec,zVal,zWrk,iFail)
C*********************************************************************
C "VecVal": search for eigen VECtors and VALues
C ..........................................................
C          .................................................
C               Part of a Package for Diaginalization of
C                        Complex Operators
C          .................................................
C
C   The parameter "Smalst" defines the smallest number such that 1
C is different from 1+Smalst. It may be Computer dependent and
C might have to be changed.
C   The parameter MxIter defines the max number of iterations
C
C Accessed Files None
C --------------
C
C Called Routines: None
C ---------------
C
C Called By: 'EigenC' 
C ---------
C
C Double Precision Complex Version
C
C J.-M. Hartmann, last change 07 March 1997
C*********************************************************************
C
      Implicit None
         Integer nTrue,nDim,K,L,iFail    
         Integer MxIter,nIter
         Integer I,J,iEnd,II,IP1,IM1,M,NN,NPL,KK,NNM1,NNM2,MM1
         Integer nnmj,nm,MM,MMM1,MP1,nm1,JJ,jm1
      Double Precision X1,X2,Y1,Y2,ZZR,YR,XR,YI,XI
      Double Precision Zero,One,Two,Norm,Smalst
         Double Complex zVal,zMat,zVec,zWrk
         Double Complex zZero,zOne,zT,zS,zX,zY,zC,zDum
      Dimension zMat(nDim,nDim),zVec(nDim,nDim)
      Dimension zVal(nDim),zWrk(nDim)
      Intrinsic dcmplx, dimag, dconjg
C
      Zero=0.d0
      One=1.d0
      Two=2.d0
      zZero=(0.d0,0.d0)
      zOne=(1.d0,0.d0)
C Smalst : the smallest (computer dependent) value such that
C 1.+Smalst is different from 1
      Smalst=1.d-14
C Max number of iterations to search Eigemvalues and vectors
      MxIter=50
C
      iFail=0
      zT=zZero
      Do 1 I=1,nTrue
      Do 10 J=1,nTrue
      zVec(I,J)=zZero
10    Continue
      zVec(I,I)=zOne
1     Continue
C
      iEnd=L-K-1
      If( iEnd.GT.0 )Then    
      Do 2 II=1,iEnd
      I=L-II
      IP1=I+1
      IM1=I-1
      Do 20 M=IP1,L
      zVec(M,I)=zMat(M,IM1)
20    Continue
      J =zWrk(I)
        If( I.NE.J )Then    
        Do 21 M=I,L
        zVec(I,M)=zVec(J,M)
        zVec(J,M)=zZero
        zVec(J,I)=zOne
21      Continue
        End If
2     Continue
      End If
      Do 3 I=1,nTrue
      If( I.GE.K.AND.I.LE.L )GoTo 3
      zVal(I)=zMat(I,I)
 3    Continue
      NN=L
C
999   Continue
      If( NN.LT.K )GoTo 150
      nIter=0
      NNM1=NN-1
      NNM2=NN-2
      If( NN.EQ.K )goto 4
998   Continue
      NPL=NN+K
      Do 40 KK=K,NNM1
      M=NPL-KK
      MM1=M-1
      If( CDABS(zMat(M,MM1)).LE.Smalst*
     *   (CDABS(zMat(MM1,MM1))+CDABS(zMat(M,M))) )GoTo 5
40    Continue
 4    Continue
      M=K
 5    Continue
      If( M.EQ.NN )GoTo 145
      If( nIter.EQ.MxIter )Then
      iFail=1
      Return
      End If
C
      If ( MOD(nIter,10).EQ.0 )GoTo 60
      zS=zMat(NN,NN)
      zX=zMat(NNM1,NN)*zMat(NN,NNM1)
      If( zX.EQ.zZero )GoTo 65
      zY=(zMat(NNM1,NNM1)-zS)/Two
      zC=CDSQRT(zY*zY+zX)
      If( DREAL(zY*DCONJG(zC)).LT.Zero )zC=-zC
      zX=zX/(zY+zC)
      zS=zS-zX
      GoTo 65
60    Continue
      X1=DREAL(zMat(NN,NNM1))
      if(nnm2.gt.0) X2=DREAL(zMat(NNM1,NNM2))  ! GCT 2012-02-18
      Y1=DIMAG(zMat(NN,NNM1))
      if(nnm2.gt.0) Y2=DIMAG(zMat(NNM1,NNM2))  ! GCT 2012-02-18
      zS=dcmplx(DABS(X1)+DABS(X2),DABS(Y1)+DABS(Y2))
65    Continue
      Do 70 I=K,NN
      zMat(I,I)=zMat(I,I)-zS
70    Continue
      zT=zT+zS
      nIter=nIter+1
C
      XR=CDABS(zMat(NNM1,NNM1))
      YR=CDABS(zMat(NN,NNM1))
      ZZR=CDABS(zMat(NN,NN))
      NNMJ=NNM1-M
      If( NNMJ.EQ.0)GoTo 80
C
      Do 75 NM=1,NNMJ
      MM=NN-NM
      MMM1=MM-1
      YI=YR
      YR=CDABS(zMat(MM,MMM1))
      XI=ZZR
      ZZR=XR
      XR=CDABS(zMat(MMM1,MMM1))
      If( YR.LE.(Smalst*ZZR/YI)*(ZZR+XR+XI) )GoTo 85
75    Continue
80    Continue
      MM=M
C
85    Continue
      MP1=MM+1
      Do 110 I=MP1,NN
      IM1=I-1
      zX=zMat(IM1,IM1)
      zY=zMat(I,IM1)
      If( CDABS(zX).GE.CDABS(zY) )GoTo 95
C
      Do 90 J=IM1,nTrue
      zDum=zMat(IM1,J)
      zMat(IM1,J)=zMat(I,J)
      zMat(I,J)=zDum
90    Continue
      zC=zX/zY
      zVal(I)=One
      GoTo 100
95    Continue
      zC=zY/zX
      zVal(I)=-One
100   Continue
      zMat(I,IM1)=zC
      Do 105 J=I,nTrue
      zMat(I,J)=zMat(I,J)-zC*zMat(IM1,J)
105   Continue
110   Continue
C                     
      Do 140 J=MP1,NN
      JM1=J-1
      zX=zMat(J,JM1)
      zMat(J,JM1)=zZero
C
      If( DREAL(zVal(J)).LE.Zero )GoTo 125
      Do 115 I=1,J
      zDum=zMat(I,JM1)
      zMat(I,JM1)=zMat(I,J)
      zMat(I,J)=zDum
115   Continue
      Do 120 I=K,L
      zDum=zVec(I,JM1)
      zVec(I,JM1)=zVec(I,J)
      zVec(I,J)=zDum
120   Continue
C
125   Continue
      Do 130 I=1,J
      zMat(I,JM1)=zMat(I,JM1)+zX*zMat(I,J)
130   Continue
      Do 135 I=K,L
      zVec(I,JM1)=zVec(I,JM1)+zX*zVec(I,J)
135   Continue
C
140   Continue
      GoTo 998
C
145   Continue
      zVal(NN)=zMat(NN,NN)+zT
      NN=NNM1
      GoTo 999
C
150   Continue
      If( nTrue.EQ.1 )Return
      Norm=Zero
      Do 160 I=1,nTrue
      Norm=Norm+CDABS(zVal(I))
      If( I.EQ.nTrue )GoTo 160
      IP1=I+1
      Do 155 J=IP1,nTrue
      Norm=Norm+CDABS(zMat(I,J))
155   Continue
160   Continue
      If( Norm.EQ.Zero )Return
C                                
      Do 180 NM=2,nTrue
      NN=nTrue+2-NM
      zX=zVal(NN)
      NNM1=NN-1
C
      Do 175 II=1,NNM1
      I=NN-II
      zC=zMat(I,NN)
      If( I.EQ.NNM1)GoTo 170
      IP1=I+1
      Do 165 J=IP1,NNM1
      zC=zC+zMat(I,J)*zMat(J,NN)
165   Continue
170   Continue
      zY=zX-zVal(I)
      If( zY.EQ.zZero )zY=Smalst*Norm
      zC=zC/zY
      zMat(I,NN)=zC
175   Continue
180   Continue
      NM1=nTrue-1
C
      Do 190 I=1,NM1
      If( I.GE.K.AND.I.LE.L )GoTo 190
      IP1=I+1
      Do 185 J=IP1,nTrue
      zVec(I,J)=zMat(I,J)
185   Continue
190   Continue
      If( L.EQ.0 )Return
      NPL=nTrue+K
      Do 200 JJ=K,NM1
      J=NPL-JJ
      Do 200 I=K,L
      zC=zVec(I,J)
      MM=J-1
      If( L.LT.J )MM=L
      Do 195 M=K,MM
      zC=zC+zVec(I,M)*zMat(M,J)
195   Continue
      zVec(I,J)=zC
200   Continue
      Return
      End
      Subroutine InvMat(zMat,zMatM1,nDim,nTrue,iFail)
C*********************************************************************
C "InvMat": INVert MATrix  
C ..........................................................
C          .   Subroutine compute the inverse of a complex .
C          .   Matrix by Gauss triangulation pivoting      .
C          .................................................
C
C Input/Output Parameters of Routine                       
C ---------------------------------
C           zMat : Complex matrix to invert (Input).
C         zMatM1 : Complex matrix=inverse of zMat (Output).
C           nDim : Number of rows (and columns of the matrix (Input)
C          nTrue : True dimension of the matrix (order) (Input)    
C          iFail : should be zero on output. =1 means a singular
C                  matrix
C
C Accessed Files None
C --------------
C
C Called Routines: None
C ---------------
C
C Called By: 'EqvLines' (eventually)
C ---------
C
C Double Precision Complex Version
C
C J.-M. Hartmann, last change 06 March 1997
C*********************************************************************
C
      Implicit None
      Integer iPivt,IndRow,IndCol      
      Integer nDim,nTrue,nLmx,iFail
      Integer I,J,K,L,LL,iRow,iCol
      Double Precision Big
      Double Complex zMat,zMatM1,zDum,zPivM1,zOne,zZero
      Parameter (nLmx=102)
      Dimension zMat(nDim,nDim),zMatM1(nDim,nDim),iPivt(nLmx)
      Dimension IndRow(nLmx),IndCol(nLmx)
      Intrinsic dcmplx
C      
C Check that dimensions are fine. Stop if not
      If ( nTrue .GT. nLmx ) Then
      Write(*,1000)nDim,nTrue       
      Stop
      End If
 1000 Format(//,1x,'************ PROBLEM !!!! ******************',
     /       /,1x,'Arrays in Routine "InvMat" are too small',
     /       /,1x,'change the value of nLmx (now ',I3,') to',I4)
C
      zOne=dcmplx(1.d0,0.D0)
      zZero=dcmplx(0.d0,0.D0)
C      
      Do 1 I=1,nTrue
      Do 1 J=1,nTrue
      zMatM1(I,J)=zMat(I,J)
1     Continue
        Do 2 J=1,nTrue
        iPivt(J)=0
2       Continue
      Do 3 I=1,nTrue
        Big=0.d0
        Do 30 J=1,nTrue
          If( iPivt(J).NE.1 )Then
            Do 300 K=1,nTrue
              If( iPivt(K).EQ.0 )Then
                If( CDABS(zMatM1(J,K)).GE.Big )Then
                  Big=CDABS(zMatM1(J,K))
                  iRow=J
                  iCol=K
                End If
              Else
              If( iPivt(K).GT.1 )Then
                iFail=1
                Return
              End If
              End If
300         Continue
          End If
30      Continue
        iPivt(iCol)=iPivt(iCol)+1
        If( iRow.NE.iCol )Then
          Do 31 L=1,nTrue
            zDum=zMatM1(iRow,L)
            zMatM1(iRow,L)=zMatM1(iCol,L)
            zMatM1(iCol,L)=zDum
31        Continue
        End If
        IndRow(I)=iRow
        IndCol(I)=iCol
           If( CDABS(zMatM1(iCol,iCol)).EQ.0.D0 )Then
           iFail=1
           Return
           End If
        zPivM1=zOne/zMatM1(iCol,iCol)
        zMatM1(iCol,iCol)=zOne
        Do 32 L=1,nTrue
          zMatM1(iCol,L)=zMatM1(iCol,L)*zPivM1
32      Continue
        Do 33 LL=1,nTrue
          If( LL.NE.iCol )Then
            zDum=zMatM1(LL,iCol)
            zMatM1(LL,iCol)=zZero
            Do 330 L=1,nTrue
              zMatM1(LL,L)=zMatM1(LL,L)-zMatM1(iCol,L)*zDum
330         Continue
          End If
33      Continue
3     Continue
      Do 4 L=nTrue,1,-1
        If( IndRow(L).NE.IndCol(L) )Then
          Do 40 K=1,nTrue
            zDum=zMatM1(K,IndRow(L))
            zMatM1(K,IndRow(L))=zMatM1(K,IndCol(L))
            zMatM1(K,IndCol(L))=zDum
40        Continue
        End If
 4    Continue
      iFail=0
      Return
      End
