      Subroutine CompAbs_qSDVFullLM(path_to_file,fname,nlev,T,P,
     &       vmr,d,fzero,grid,ncp,nspeci,ispeci,mm,molwet,nu1,nu2,
     &       vvoigt,vac)

C    Define Variables
      integer*4 ncp,nlev,ispeci,nspeci,mm,ns,mp,np,iFile,lp,
     &molwet,i,j,kv1,kv2,nv,n,jv,Wfile,
     &Wfile2,WTfile,WTfile2,Fail
      parameter (mp=50)
      real*4 T(nlev),P(nlev),vmr(nspeci,nlev),d(nlev),vvoigt(ncp),
     &vac(ncp,nlev,0:*),con,y,T0,dopphwm
      real*8 fzero,grid,nu1,nu2,abhw(mp),Tabhw(mp),apSHFT(mp),
     &vo(mp),E(mp),sbhw(mp),spSHFT(mp),S(mp),SD(mp),SD_s(mp),
     &dn(mp),TapSHFT(mp),TspSHFT(mp),dgen(mp),
     &Ro(mp),Ro0(mp),QT,QT0,Trat,Td,frac,pi,srpi,tnulst,vcent,fmax,godw,
     &sg,sxcgsoripdw,flinwid,mat(mp,8),W_air(mp,mp),WT_air(mp,mp),
     &W_self(mp,mp),WT_self(mp,mp),x,GG(mp,4),temp(1,4)
      complex*16 W(mp,mp),H(mp),A(mp,mp),AT(mp,mp),Wrk(mp),z1,z2
     &,G(mp)
      character path_to_file*(*),WaFile*10,WsFile*11,WaTFile*11,
     &WsTFile*12,cnp*2,fname*14
      parameter (T0=296.0,pi=3.141593d0,con=1.43881,
     &           Wfile=1,WTfile=2,iFile=3,Wfile2=4,WTfile2=5)

C    Calculate constants and set variables outside loop over levels
      srpi=dsqrt(pi)
      frac=5.0E-06
      tnulst=1.0E-13
      fmax=fzero+ncp*grid
      !calculate the partion function at 296 K
      QT0=Qmol(T0,ispeci)

C    Read the first two line of the File with the spec parameters
      lp=len_trim(path_to_file)
      open(Unit=iFile,File=path_to_file(:lp)//fname
     &,Status='old')
      Read(iFile,'(I3,1x,I3)',End=100)ns,np
      if(np.gt.mp) then
         write(*,*) ' np, mp = ',np,mp
         stop 'Cal_FullLM: Increase parameter mp '
      endif
      Read(iFile,'(a10,1x,a11,1x,a11,1x,a12)',End=100)WaFile,WaTFile,
     &WsFile,WsTFile
C    Loop to read the specparameters from the file
      Do i=1,np
         Read(iFile,1001,End=100)vo(i),S(i),abhw(i),Tabhw(i),E(i),
     &   apSHFT(i),SD(i),sbhw(i),spSHFT(i),SD_s(i),dn(i),
     &   TapSHFT(i),TspSHFT(i),dgen(i)
         !Calculate the lower state population
         Ro0(i)=dgen(i)*dexp(-con*E(i)/T0)/QT0
      end do !End Loop over reading lines
1001   Format(f12.6,E10.3,f6.4,f5.3,f10.4,f9.6,f6.3,f7.5,f8.5,
     &        f8.5,f8.5,f9.6,f9.6,88x,f6.1)
100   Close(iFile)

C     Open and read in W matrix elements into W matrix
      open(Unit=Wfile,File=path_to_file(:lp)//WaFile,Status='old')
      open(Unit=WTfile,File=path_to_file(:lp)//WaTFile,Status='old')
      open(Unit=Wfile2,File=path_to_file(:lp)//WsFile,Status='old')
      open(Unit=WTfile2,File=path_to_file(:lp)//WsTFile,Status='old')
      write(cnp,'(I2)')np !Convert np to char
c      write(*,*) 'np=',np,'('//cnp//'f9.6)'
      Do i=1,np 
c         read(Wfile,'('//cnp//'f9.6)')W_air(i,:)
c         read(WTfile,'('//cnp//'f5.2)')WT_air(i,:)
c         read(Wfile2,'('//cnp//'f9.6)')W_self(i,:)
c         read(WTfile2,'('//cnp//'f5.2)')WT_self(i,:)
         read(Wfile,'('//cnp//'f9.6)')(W_air(i,j),j=1,np)
         read(WTfile,'('//cnp//'f5.2)')(WT_air(i,j),j=1,np)
         read(Wfile2,'('//cnp//'f9.6)')(W_self(i,j),j=1,np)
         read(WTfile2,'('//cnp//'f5.2)')(WT_self(i,j),j=1,np)
      end do !End loop over reading W matrix elements
      close(Wfile)
      close(WTfile)
      close(Wfile2)
      close(WTfile2)

C     Calculate the absorption coefficients and store in vac
      do n=1,nlev !Loop over all levels
        QT=Qmol(T(n),ispeci)
        Trat=T0/T(n)
        Td=T(n)-T0
C       Loop over line parameters
        do i=1,np
          Ro(i)=dgen(i)*dexp(-con*E(i)/T(n))/QT
          !Line Strength at given temp
          mat(i,1)=S(i)*(Ro(i)/Ro0(i))*(1-dexp(-con*vo(i)/T(n)))/
     &                                 (1-dexp(-con*vo(i)/T0))
          !Dipole Matrix elemnt at T
          mat(i,2)=dsqrt(mat(i,1)/Ro(i))
          !Lorentz width
          mat(i,3)=(abhw(i)*(1-vmr(ispeci,n))+sbhw(i)*vmr(ispeci,n))
     &             *(Trat**Tabhw(i))*P(n)
          !Pressure Shift
          mat(i,4)=((apSHFT(i)+TapSHFT(i)*Td)*(1-vmr(ispeci,n))+
     &              (spSHFT(i)+TspSHFT(i)*Td)*vmr(ispeci,n))*P(n)
          !Speed-dependent Lorentz width
          mat(i,5)=SD(i)*mat(i,3)
          !Dicke Narrowing parameter
          mat(i,6)=dn(i)*P(n)
          !Line center frequincies
          mat(i,8)=vo(i)
        end do
        if(P(n).gt.0.0)then
C       Create the Relaxtion matrix
        A(:,:)=dcmplx(0.d0,0.d0)
        AT(:,:)=dcmplx(0.d0,0.d0)
        H(:)=dcmplx(0.d0,0.d0)
        Wrk(:)=dcmplx(0.d0,0.d0)
        W(:,:)=dcmplx(0.d0,0.d0)
        do i=1,np
           do j=1,np
              if(W_air(i,j).gt.0.0d0)then
                 x=-P(n)*(W_air(i,j)*(trat**WT_air(i,j))*(1
     &             -vmr(ispeci,n))+W_self(i,j)*(trat**WT_self(i,j))
     &             *(vmr(ispeci,n)))
                 W(i,j)=dcmplx(0.0d0,x)
                 W(j,i)=dcmplx(0.0d0,imag(W(i,j))*Ro(j)/Ro(i))
              end if
           end do
           W(i,i)=dcmplx(vo(i),mat(i,3))
        end do
C       Calculate the Effective Line parameters using the W matrix
        call EigenC(np,mp,W,H,A,Wrk,Fail)
        call InvMat(A,AT,mp,np,Fail)
        do i=1,np
           z1=0.d0
           z2=0.d0
           do j=1,np
              z1=z1+Ro(j)*mat(j,2)*A(j,i)
              z2=z2+mat(j,2)*AT(i,j)
           end do
           G(i)=z1*z2
           GG(i,1)=real(G(i))
           GG(i,2)=imag(G(i))/real(G(i))
           GG(i,3)=real(H(i))
           GG(i,4)=imag(H(i))
         !write(*,*)real(G(i))
        end do
C       Sort the lines so they match the SD width and Shift calculated
        do i=1,np-1
           do j=i+1,np
              if (GG(i,3).gt.GG(j,3))then
                  temp(1,:)=GG(j,:)
                  GG(j,:)=GG(i,:)
                  GG(i,:)=temp(1,:)
              end if
           end do
        end do
C       Replace the line parameters with the parameters calculated using
C       Full Line Mixing
        do i=1,np
           mat(i,1)=GG(i,1)
           mat(i,3)=GG(i,4)
           mat(i,7)=GG(i,2)
           if(abs(GG(i,2)).le.0.0d0)then
              mat(i,8)=mat(i,8)
           else
              mat(i,8)=GG(i,3)
           endif
        !write(*,*)i,mat(i,8),mat(i,1)
        end do
        end if !end if P(n)>0

C       Loop over line parameters to Calculate VAC
        do i=1,np
C          Skip lines outside of Transmission calculation
           if(mat(i,8).lt.nu1)goto 40
           if(mat(i,8).gt.nu2)goto 44
           !Doppler width
           dopphwm=4.3014e-7*sngl(mat(i,8))*sqrt(T(n)/molwet)
           !Calculate the start and last index of the grid
           sxcgsoripdw=mat(i,1)*d(n)*vmr(ispeci,n)*(1/srpi)/dopphwm
           godw=grid/dopphwm
           y=sngl(mat(i,3)/dopphwm)
           flinwid=dble((3+amin1(y/frac,sqrt(abs(y*sxcgsoripdw/srpi
     &                                            /tnulst))))/godw)
           vcent=(mat(i,8)-fzero+dble(mat(i,4)))/grid
           !Start index
           kv1=max0(1,1+int(vcent-flinwid))
           !End index
           kv2=min0(int(vcent+flinwid),ncp)
           !Number of grid points
           nv=kv2-kv1+1
           !Calculate the line shape for each point
           if(mat(i,8).ge.fmax) then
              mat(:,:)=0.d0
              goto 44
           else if((kv1.ge.1).and.(kv2.le.ncp))then
              sg=(dble(kv1)*grid)+fzero
              call qSDV(mat(i,8),dopphwm,mat(i,3),mat(i,5),mat(i,4),
     &                  0.0d0,sg,mat(i,7),nv,grid,vvoigt(kv1))
              !Put abs into vac
              do jv=kv1,kv2
                 vac(jv,n,mm)=vac(jv,n,mm)-sngl(sxcgsoripdw)*vvoigt(jv)
              end do
           end if
  40       continue
        end do !End Loop to calculate VAC
  44   continue
      end do !End loop over levels

      End !End subroutine

C******************************************************

      REAL function Qch4(iso,T)
c ch4 partition function
      implicit none
      Integer*4 iso,j
      Real*4 T
      Real*8 ch4_Qcoef(3,4)
c
c...ch4  --  211
      DATA (ch4_Qcoef(1,j),j=1,4)/-.17475d+02, .95375d+00,
     +                .39758d-02,-.81837d-06/
c...ch4  --  311
      DATA (ch4_Qcoef(2,j),j=1,4)/-.27757d+02, .17264d+01,
     +                .93304d-02,-.48181d-05/
c...ch4  --  212
      DATA (ch4_Qcoef(3,j),j=1,4)/-.89810d+03, .44451d+02,
     +                .17474d+00,-.22469d-04/
      if( (T.lt.70.d0) .or. (T.gt.415.d0) )then
      write(6,*)T,' out of range for ch4 partition function'
      stop
        else
        Qch4 = sngl(ch4_Qcoef(iso,1)
     +       +      ch4_Qcoef(iso,2)*T
     +       +      ch4_Qcoef(iso,3)*T*T
     +       +      ch4_Qcoef(iso,4)*T*T*T)
        !write(*,*)Qch4
      endif
      RETURN
      END
C------------------------------------------------------
C******************************************************

      REAL function Qco2(T)
C CO2 partition function where polynomial values determined
C by fitting Q of CO2 (main iso), where Q was calculated
C using TIPS 2008 for T = 70-415 K.  
      implicit none
      Real*4 T
c
      if( (T.lt.70.d0) .or. (T.gt.415.d0) )then
        write(6,*)T,' out of range for ch4 partition function'
        stop
      else
        Qco2 = -1.351258
     +       + (0.95261)*T
     +       - (0.0007252)*T*T
     +       + (2.6943e-6)*T*T*T
      endif
      RETURN
      END
C------------------------------------------------------
C******************************************************

      REAL function Qmol(T,ispec)
C CO2 partition function where polynomial values determined
C by fitting Q of CO2 (main iso), where Q was calculated
C using TIPS 2008 for T = 70-415 K.
      implicit none
      Integer*4 ispec
      Real*4 T
c
      Qmol=0.0
      if( (T.lt.70.d0) .or. (T.gt.415.d0) )then
        write(6,*)T,' out of range for ch4 partition function'
        stop
      else
c       CO2 mol=2 iso=1      
        if(ispec.eq.4)Qmol = -1.351258+(0.95261)*T-(0.0007252)*T**2
     +                     +(2.6943e-6)*T**3
C       CH4 mol=6 iso=1 DW changed to 35 to match new isotopologs file
        if(ispec.eq.35)Qmol = -.17475e+02+.95375e+00*T
     +                      +(.39758e-02)*T**2+(-.81837e-06)*T**3
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
      !Zvec(:,:)=(0.0d0,0.0d0)
      !zVal(:)=(0.0d0,0.0d0)
      !zWrk(:)=(0.0d0,0.0d0)
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
        Norm=Norm+dreal(zVec(iRow,iCol)*DCONJG(zVec(iRow,iCol)))
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
      KK=int(DIMAG(zWrk(I)))
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
      If( CDABS(zMat(I,J)).gt.Zero )GoTo 40
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
      If( CDABS(zMat(I,J)).gt.Zero )GoTo 60
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
      If( abs(zX-zZero).le.0.0 )GoTo 1
      MP1=M+1
      Do 13 I=MP1,L
      zY=zMat(I,M-1)
      If( abs(zY-zZero).gt.0.0 )Then
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
      zS=(0.0,0.0) ! Avoid compiler warning (may be used unitialized)
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
      J =int(dreal(zWrk(I)))
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
      If( abs(zX-zZero).le.0.0 )GoTo 65
      zY=(zMat(NNM1,NNM1)-zS)/Two
      zC=CDSQRT(zY*zY+zX)
      If( DREAL(zY*DCONJG(zC)).LT.Zero )zC=-zC
      zX=zX/(zY+zC)
      zS=zS-zX
      GoTo 65
60    Continue
      if(nnm2.gt.0) then
         X1=DREAL(zMat(NN,NNM1))
         X2=DREAL(zMat(NNM1,NNM2))
         Y1=DIMAG(zMat(NN,NNM1))
         Y2=DIMAG(zMat(NNM1,NNM2))
         zS=dcmplx(DABS(X1)+DABS(X2),DABS(Y1)+DABS(Y2))
      end if
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
      If( abs(Norm-Zero).le.0.0 )Return
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
      If( abs(zY-zZero).le.0.0 )zY=Smalst*Norm
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

      iRow=0 ! Do-nothing statement to prevent compiler warning (may be used uninitialized)
      iCol=0 ! Do-nothing statement to prevent compiler warning (may be used uninitialized)
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
           If( CDABS(zMatM1(iCol,iCol)).le.0.0D0 )Then
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
