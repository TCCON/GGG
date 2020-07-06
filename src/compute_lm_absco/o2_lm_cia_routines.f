C*********************************************************************
      Subroutine LMandCIAO2(path_to_input_files,T,Ptot,
     &  xo2,xh2o,SigMin,DSig,nSig, AbsV, AbsY, AbsW)
C*********************************************************************
C"LMandCIAO2": COMPute ABSorptions
C.................................................
C.     Subroutine to Compute the absorption      .
C.     Coefficient due to all Lines retained     .
C.     by Subroutine "Readline"                  . 
C.     Subroutine "Readline".                    .
C.     Voigt, First Order L-M                    .
C.................................................
C
C   Inputs:
C       path_to_input_files: where to find linelist, matrices, etc.
C    T       : Temperature in Kelvin 
C    Ptot    : Total pressure in Atmosphere 
C    xO2     : O2 volume mixing ratio 
C    xH2O    : H2O vmr
C    SigMin  : Minimum WaveNumber of the Computation (Cm-1)
CDSig    : WaveNumber Step of the Computation (Cm-1)
C       nSig    : Number of frequencies
C
C
C   Outputs
CAbsV    : Absorption Coefficient neglecting LineMixing
C                (assuming Voigt Line-Shapes) (cm-1)
CAbsY    : Difference between Absorption Coefficient predicted
C                using First Order Line-Mixing and Voigt (cm-1)
C       AbsW    : Absorption Coefficient due to CIA (cm-1)
C
C
C   Other important Input Quantities (through Common Statements)
C--------------------------------
CHWT     : Air Broadened HalfWidths of the Lines for the
C          Considered Temperature and Pressure (Cm-1)
Co2_PopuT   : Populations of the Lower Levels of the Lines
C  at Temperature T
Co2_YT      : Air Broadened First Order Line Mixing Coefficients 
C               of the Lines for the Considered Temperature and
C          Pressure (No Unit)
CDipo    : Dipole transition Moments of the Lines
C
CAccessed Files:  None
C--------------
C
CCalled Routines: 'ConvTP' (CONVert data to current T and Press)
C---------------  'CPF'  (Complex Probability Function)
C
CCalled By: Main Program
C---------
C
C     Double Precision Version
C
C     H. Tran, last change 24 November 2005
C*********************************************************************
C
      Implicit None
      include '../gfit/ggg_int_params.f'
      Integer nLmx,nSigmx,nCIAmx
      Integer o2_nLines,nCIA
      Integer nSig,iSig,iLine
      Integer ind,idum
      Double Precision T,Ptot,xO2,xh2o,SigMin,DSig
      Double Precision o2_Sig,Dipo,HWT,o2_PopuT,o2_YT
      Double Precision aMass,Pi
      Double Precision SigMoy,GamD,Cte,Cte1,Cte2
      Double Precision SigC
      Double Precision aa,bb
      Double Precision XX,YY,WR,WI
      Double Precision SigCIA,CIA,Acialoc
      Double complex zint,zval
      character path_to_input_files*(mpath+80)
C Max Number of Lines, of Spectral and of CIA data points
      Parameter (nLmx=100)
      Parameter (nSigmx=500000)
      Parameter (nCIAmx=100)
C Characteristic of the Band
      Double Precision AbsV(nsigmx),AbsY(nsigmx),AbsW(nsigmx)
      Common/o2_Bands/o2_nLines,nCIA
      Common/Mean/SigMoy
C Data of the Lines
      Common/o2_LineSg/o2_Sig(nLmx) 
      Common/o2_GamT/HWT(nLmx) 
      Common/DipolT/Dipo(nLmx) 
      Common/o2_PopuT/o2_PopuT(nLmx)
      Common/o2_YLT/o2_YT(nLmx)
      Common/o2_eqvLines/zval(nlmx),zint(nlmx)
C Data of CIA part
      Common/o2_LineSgCIA/SigCIA(nCIAmx)
      Common/SCIA/CIA(nCIAmx)
C Results (Absorption Coefficients)
c      Common/CabsVO2/AbsV(nSigmx)
c      Common/CabsFO2/AbsY(nSigmx)
C Constants      
      Data Pi/3.141592654d0/
      Data aMass/32.D-3/

      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mauxcol ! Avoid compiler warning (unused parameter)
      idum=mcolvav ! Avoid compiler warning (unused parameter)
      idum=mgas    ! Avoid compiler warning (unused parameter)
      idum=mlev    ! Avoid compiler warning (unused parameter)
      idum=mrow_qc ! Avoid compiler warning (unused parameter)
      idum=mspeci  ! Avoid compiler warning (unused parameter)
      idum=mvmode  ! Avoid compiler warning (unused parameter)
      idum=ncell   ! Avoid compiler warning (unused parameter)
      idum=nchar   ! Avoid compiler warning (unused parameter)

C----------
C
      xx=xh2o  ! avoid compiler warning (unused variable)
C Convert Band Data to Considered Temperature and Pressure
      Call ConvTPO2(path_to_input_files,T,Ptot)
C
C Compute the Doppler Width at SigMoy
      GamD=DSQRT(2.d0*DLOG(2.d0)*8.314d0/(3d+8**2.0))*
     1 DSQRT(T/aMass)*SigMoy
      Cte=DSQRT(DLOG(2.d0))/GamD
      Cte1=Cte/DSQRT(Pi)
c set to zero
      do iSig=1,nSig
         AbsV(isig)=0.D0
         AbsY(isig)=0.D0
         AbsW(isig)=0.D0
      end do
C
C Various Lines
C
      do 3 iLine=1,o2_nLines
      SigC=SigMin-DSig
C
C Various WaveNumber
      Do 2 iSig=1,nSig
      SigC=SigC+DSig
C
      Cte2=SigC*(1.d0-DEXP(-1.4388*SigC/T))
C Complex Probability Function for the "Real" Q-Lines
      XX=(o2_Sig(iLine)-SigC)*Cte
      YY=HWT(iLine)*Ptot*Cte
      Call CPF(XX,YY,WR,WI)
C Voigt Absorption Coefficient
      aa=Cte1*xO2*Ptot*o2_PopuT(iLine)*(Dipo(iLine)**2)*WR
      AbsV(isig)=AbsV(isig)+aa*Cte2
C First Order Line-Mixing Absorption Coefficient
c      bb=Cte1*xO2*Ptot*o2_PopuT(iLine)*(Dipo(iLine)**2)*(WR-o2_YT(iLine)*WI)
      bb=-Cte1*xO2*Ptot*o2_PopuT(iLine)*(Dipo(iLine)**2)*WI*o2_YT(iLine)
      AbsY(isig)=AbsY(isig)+bb*Cte2
2     Continue
3     Continue
C 
C Calculation of the absorption coefficient due to CIA contribution
C
      SigC=Sigmin-Dsig
      do iSig=1,nSig
      SigC=SigC+dSig
      ind=(SigC-SigCIA(1))/10.+1
      if(ind.le.0 .or. ind.gt.nCIAmx) then
         write(*,*)'Warning: Outside region for O2 CIA'
         Acialoc=0.0
      else
         Acialoc=CIA(ind)+(SigC-SigCIA(ind))
     &    *(CIA(ind+1)-CIA(ind))/10.
      endif
        AbsW(isig)=Acialoc*xo2
c       AbsY(isig)=AbsY(isig)+Acialoc*xo2
      enddo

      Return
      End Subroutine LMandCIAO2
C
C -------------------------------------------------
C*********************************************************************
      Subroutine ConvTPO2(path_to_input_files,T,Ptot)
C*********************************************************************
C      "ConvTP": CONVert to Temperature and Pressure
C........................................................
C      .Subroutine to convert the data read by Subroutine  .
C      ."Readline" to the condition of temperature 'T' and .
C      .pressure 'Ptot' for all retained lines             . 
C........................................................
C
C      Input/Output Parameters of Routine (Arguments or Common)
C----------------------------------
C      o2_nLines  : Number of lines of A band (Input).
C      T    : Temperature in Kelvin (Input).
C      Ptot    : Total Pressure in Atmosphere (Input).
C      SigMoy  : The Population-Averaged value of the Positions
C               of the Lines in the current band (Output).
C
C      Other important Output Quantities (through Common Statements)
C----------------------------------
C      HWT: Air Broadened HalfWidths of the Lines for the
C           Considered Temperature and Pressure (Cm-1)
C      o2_PopuT   : Populations of the Lower Levels of the Lines
C           at Temperature T
C      o2_YT      : Air Broadened First Order Line Mixing Coefficients 
C           of the Lines for the Considered Temperature and
C           Pressure (No Unit)
C ..........................................................
C
C     Accessed Files:  None
C     --------------
C     
C     Called Routines: "PFO2" (Partition Function of O2)
C     ---------------  "Readline" (READ LINE data and relaxation matrix)
C     
C     Called By: "CompAbs"(COMPute ABSorpton)
C     ---------
C     
C     Double Precision Version
C     
C     H. Tran, last change 28 November 2005
C*********************************************************************
C
      Implicit double precision (a-h,o-z)
      include '../gfit/ggg_int_params.f'
      Integer o2_nLines,idum
      Double complex zInt,zVal
      character path_to_input_files*(mpath+80)
C Max Number of Lines and of matrix elements
      Parameter (nLmx=100)
      Parameter (nMmx=5000)
      Parameter (nCIAmx=100)
C Constants
      Common/Cst/T0,A
C Characteristic of the Band
      Common/o2_Bands/o2_nLines,nCIA
      Common/Mean/SigMoy
C Data of CIA part
      Common/o2_LineSgCIA/SigCIA(nCIAmx)
      Common/SCIA/CIA(nCIAmx)
C Lines data at Ref Temperature/Pressure
      Common/o2_LineSg/o2_Sig(nLmx)
      Common/o2_GamT0/HWT0(nLmx) 
      Common/DTGAM/BHW(nLmx) 
      Common/o2_PopTrf/o2_PopuT0(nLmx)
      Common/o2_Energy/E(nLmx) 
      Common/DipolT0/Dipo0(nLmx) 
      Common/DipolR/DipoR(nLmx)
      Common/RelxT0/WT0(nMmx)
      Common/dTRelx/BW(nMmx)
      Common/DEL/SHIFT(nLmx)
      Common/NJline/m1(nLmx)
C Lines data at (Temperature,Pressure) 
      Common/o2_GamT/HWT(nLmx) 
      Common/o2_PopuT/o2_PopuT(nLmx)
      Common/DipolT/Dipo(nLmx) 
      Common/o2_YLT/o2_YT(nLmx)
      Common/o2_eqvLines/zval(nlmx),zint(nlmx)
C Other parameters
      Double precision WW0(nLmx,nLmx)
      Double Precision Sup(nLmx)
      Double Precision Slow(nLmx)
c      DImension zVec(nlmx,nlmx),zVecM1(nlmx,nlmx)
C----------
      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mauxcol ! Avoid compiler warning (unused parameter)
      idum=mcolvav ! Avoid compiler warning (unused parameter)
      idum=mgas    ! Avoid compiler warning (unused parameter)
      idum=mlev    ! Avoid compiler warning (unused parameter)
      idum=mrow_qc ! Avoid compiler warning (unused parameter)
      idum=mspeci  ! Avoid compiler warning (unused parameter)
      idum=mvmode  ! Avoid compiler warning (unused parameter)
      idum=ncell   ! Avoid compiler warning (unused parameter)
      idum=nchar   ! Avoid compiler warning (unused parameter)

      T0=296.D0
      A=1.4388D0
      Amagat=Ptot*273.15/T
      call QT_O2(T,1,gsi,QT)                      
      call QT_O2(T0,1,gsi,QT0)
      rapq1=qt0/QT
C Call "Readline" subroutine
      Call ReadDataO2(path_to_input_files)
C Collisional shift 
      do iLine=1,o2_nLines
      o2_Sig(iLine)=o2_Sig(iLine)+Ptot*SHIFT(iLine)
      enddo
C Calculation of the population of the lower level, dipole elements 
C and relaxation matrix elements at temperature T
      do iLine=1,o2_nLines
      o2_PopuT(iLine)=o2_PopuT0(iLine)*rapq1
     & *DEXP(A*E(iLine)*(1/T0-1/T))
      Dipo(iLine)=Dipo0(iLine)*DSQRT(T0/T)
      do iLineP=1,o2_nLines
      LK=iLineP+(iLine-1)*o2_nLines
      WW0(iLine,iLineP)=WT0(LK)*(T0/T)**BW(LK)
      enddo
      enddo
C Calculation of matrix elements of the upward transitions from the detailed balance
      do 225 iLine=1,o2_nLines
      do 227 iLineP=1,o2_nLines
      if (E(iLineP).lt.E(iLine)) goto 227
          WW0(iLineP,iLine)=WW0(iLine,iLineP)*o2_PopuT(iLineP)/
     &    o2_PopuT(iLine)
227      continue
225      continue
C Diagonal elements by line-broadening at temperature T
      do iLine=1,o2_nLines
      HWT(iLine)=HWT0(iLine)*(T0/T)**BHW(iLine)
      WW0(iLine,iLine)=HWT(iLine)
      enddo
C Sigmoy
      SigMoy=0.d0
      SumWgt=0.d0
      do iLine=1,o2_nLines
      Wgt=o2_PopuT(iLine)*(Dipo(iLine)**2)
      SigMoy=SigMoy+o2_Sig(iLine)*Wgt
      SumWgt=SumWgt+Wgt
      enddo
      SigMoy=SigMoy/SumWgt
C****************************************************
C "Renormalization" procedure
C
      Do 295 iLine=1,o2_nLines
      Sup(iLine)=0.d0
      Slow(iLine)=0.d0
      Do 297 iLineP=1,o2_nLines
      if (iLineP.le.iLine) then
      Sup(iLine)=Sup(iLine)+DipoR(iLineP)*WW0(iLineP,iLine)
      else
      Slow(iLine)=Slow(iLine)+DipoR(iLineP)*WW0(iLineP,iLine)
      endif
297      Continue
      Do 299 iLineP=1,o2_nLines
      if (iLineP.le.iLine) goto 299
      WW0(iLineP,iLine)=-WW0(iLineP,iLine)*
     1 (Sup(iLine)-DipoR(iLine)*(1.-abs(m1(iLine))/36.)**2.*0.04)
     2 /Slow(iLine)

        WW0(iLine,iLineP)=WW0(iLineP,iLine)*o2_PopuT(iLine)/
     &  o2_PopuT(iLineP)
299   Continue
295   Continue
C Computation of the first order line-coupling coefficients 
      do 125 iLine=1,o2_nLines
      YY=0.D0
      do iLineP=1,o2_nLines
      if (iLine.ne.iLineP.and.o2_Sig(iLine).ne.o2_Sig(iLineP)) then
       YY=YY+2.d0*Dipo(iLineP)/Dipo(iLine)*
     1 WW0(iLineP,iLine)/(o2_Sig(iLine)-o2_Sig(iLineP))
      endif
      enddo
      o2_YT(iLine)=Ptot*YY


125   continue
C
C Conversion of the CIA data
C
      do iCIA=1,nCIA
      CIA(iCIA)=CIA(iCIA)*Amagat**2
      enddo
      Return
      End Subroutine ConvTPO2

C
C -------------------------------------------------
C*********************************************************************
      Subroutine QT_O2    (                       
     I T,     ! temperature in K 
     I iso,   ! isotope code (HITRAN INDEX)
     O gsi,   ! state independent nuclear degeneracyfactor
     O QT)    ! Total Internal Partition Function
 
      implicit DOUBLE PRECISION (a-h,o-z)
      parameter (NT=119)
      
      dimension xgj( 3), QofT( 3,119),Q(NT),tdat(NT)
      data xgj/ 1.,1.,6./
      data Tdat/  60.,  85., 110., 135., 160., 185., 210., 235.,
     + 260., 285., 310., 335., 360., 385., 410., 435., 460., 485.,
     + 510., 535., 560., 585., 610., 635., 660., 685., 710., 735.,
     + 760., 785., 810., 835., 860., 885., 910., 935., 960., 985.,
     +1010.,1035.,1060.,1085.,1110.,1135.,1160.,1185.,1210.,1235.,
     +1260.,1285.,1310.,1335.,1360.,1385.,1410.,1435.,1460.,1485.,
     +1510.,1535.,1560.,1585.,1610.,1635.,1660.,1685.,1710.,1735.,
     +1760.,1785.,1810.,1835.,1860.,1885.,1910.,1935.,1960.,1985.,
     +2010.,2035.,2060.,2085.,2110.,2135.,2160.,2185.,2210.,2235.,
     +2260.,2285.,2310.,2335.,2360.,2385.,2410.,2435.,2460.,2485.,
     +2510.,2535.,2560.,2585.,2610.,2635.,2660.,2685.,2710.,2735.,
     +2760.,2785.,2810.,2835.,2860.,2885.,2910.,2935.,2960.,2985.,
     +3010./
c...       O2
c...        --        66
      data (QofT( 1,J),J=1,119)/ 0.44334E+02, 0.62460E+02, 0.80596E+02,
     + 0.98738E+02, 0.11688E+03, 0.13503E+03, 0.15319E+03, 0.17136E+03,
     + 0.18954E+03, 0.20775E+03, 0.22600E+03, 0.24431E+03, 0.26270E+03,
     + 0.28119E+03, 0.29981E+03, 0.31857E+03, 0.33750E+03, 0.35662E+03,
     + 0.37594E+03, 0.39550E+03, 0.41529E+03, 0.43535E+03, 0.45568E+03,
     + 0.47630E+03, 0.49722E+03, 0.51844E+03, 0.53998E+03, 0.56185E+03,
     + 0.58406E+03, 0.60660E+03, 0.62949E+03, 0.65274E+03, 0.67635E+03,
     + 0.70031E+03, 0.72465E+03, 0.74936E+03, 0.77444E+03, 0.79990E+03,
     + 0.82574E+03, 0.85197E+03, 0.87858E+03, 0.90558E+03, 0.93297E+03,
     + 0.96076E+03, 0.98895E+03, 0.10175E+04, 0.10465E+04, 0.10759E+04,
     + 0.11057E+04, 0.11359E+04, 0.11665E+04, 0.11976E+04, 0.12290E+04,
     + 0.12609E+04, 0.12931E+04, 0.13258E+04, 0.13590E+04, 0.13925E+04,
     + 0.14265E+04, 0.14609E+04, 0.14958E+04, 0.15311E+04, 0.15669E+04,
     + 0.16031E+04, 0.16397E+04, 0.16768E+04, 0.17144E+04, 0.17524E+04,
     + 0.17909E+04, 0.18298E+04, 0.18692E+04, 0.19091E+04, 0.19495E+04,
     + 0.19904E+04, 0.20318E+04, 0.20736E+04, 0.21160E+04, 0.21588E+04,
     + 0.22022E+04, 0.22461E+04, 0.22905E+04, 0.23354E+04, 0.23809E+04,
     + 0.24268E+04, 0.24734E+04, 0.25204E+04, 0.25680E+04, 0.26162E+04,
     + 0.26649E+04, 0.27142E+04, 0.27641E+04, 0.28145E+04, 0.28655E+04,
     + 0.29171E+04, 0.29693E+04, 0.30221E+04, 0.30755E+04, 0.31295E+04,
     + 0.31841E+04, 0.32393E+04, 0.32951E+04, 0.33516E+04, 0.34087E+04,
     + 0.34665E+04, 0.35249E+04, 0.35839E+04, 0.36436E+04, 0.37040E+04,
     + 0.37650E+04, 0.38267E+04, 0.38891E+04, 0.39522E+04, 0.40159E+04,
     + 0.40804E+04, 0.41455E+04, 0.42114E+04, 0.42780E+04, 0.43452E+04,
     + 0.44132E+04/
c...        --        68
      data (QofT( 2,J),J=1,119)/ 0.89206E+02, 0.12759E+03, 0.16600E+03,
     + 0.20442E+03, 0.24285E+03, 0.28128E+03, 0.31973E+03, 0.35821E+03,
     + 0.39672E+03, 0.43530E+03, 0.47398E+03, 0.51281E+03, 0.55183E+03,
     + 0.59108E+03, 0.63062E+03, 0.67051E+03, 0.71078E+03, 0.75148E+03,
     + 0.79265E+03, 0.83435E+03, 0.87659E+03, 0.91941E+03, 0.96285E+03,
     + 0.10069E+04, 0.10517E+04, 0.10971E+04, 0.11432E+04, 0.11901E+04,
     + 0.12377E+04, 0.12861E+04, 0.13352E+04, 0.13851E+04, 0.14358E+04,
     + 0.14872E+04, 0.15395E+04, 0.15926E+04, 0.16466E+04, 0.17013E+04,
     + 0.17569E+04, 0.18134E+04, 0.18706E+04, 0.19288E+04, 0.19877E+04,
     + 0.20476E+04, 0.21083E+04, 0.21698E+04, 0.22323E+04, 0.22956E+04,
     + 0.23598E+04, 0.24248E+04, 0.24908E+04, 0.25576E+04, 0.26253E+04,
     + 0.26940E+04, 0.27635E+04, 0.28339E+04, 0.29052E+04, 0.29775E+04,
     + 0.30506E+04, 0.31247E+04, 0.31997E+04, 0.32756E+04, 0.33524E+04,
     + 0.34302E+04, 0.35089E+04, 0.35885E+04, 0.36691E+04, 0.37506E+04,
     + 0.38331E+04, 0.39166E+04, 0.40010E+04, 0.40864E+04, 0.41727E+04,
     + 0.42601E+04, 0.43484E+04, 0.44377E+04, 0.45280E+04, 0.46193E+04,
     + 0.47116E+04, 0.48049E+04, 0.48992E+04, 0.49946E+04, 0.50909E+04,
     + 0.51883E+04, 0.52868E+04, 0.53863E+04, 0.54868E+04, 0.55884E+04,
     + 0.56911E+04, 0.57949E+04, 0.58997E+04, 0.60056E+04, 0.61126E+04,
     + 0.62207E+04, 0.63298E+04, 0.64401E+04, 0.65516E+04, 0.66641E+04,
     + 0.67778E+04, 0.68926E+04, 0.70085E+04, 0.71256E+04, 0.72439E+04,
     + 0.73633E+04, 0.74839E+04, 0.76056E+04, 0.77286E+04, 0.78527E+04,
     + 0.79781E+04, 0.81046E+04, 0.82324E+04, 0.83613E+04, 0.84915E+04,
     + 0.86229E+04, 0.87556E+04, 0.88895E+04, 0.90247E+04, 0.91611E+04,
     + 0.92988E+04/
c...        --        67
      data (QofT( 3,J),J=1,119)/ 0.52071E+03, 0.74484E+03, 0.96908E+03,
     + 0.11934E+04, 0.14177E+04, 0.16422E+04, 0.18667E+04, 0.20913E+04,
     + 0.23161E+04, 0.25413E+04, 0.27671E+04, 0.29936E+04, 0.32212E+04,
     + 0.34501E+04, 0.36806E+04, 0.39130E+04, 0.41476E+04, 0.43846E+04,
     + 0.46242E+04, 0.48668E+04, 0.51125E+04, 0.53615E+04, 0.56140E+04,
     + 0.58701E+04, 0.61300E+04, 0.63938E+04, 0.66617E+04, 0.69337E+04,
     + 0.72099E+04, 0.74904E+04, 0.77754E+04, 0.80647E+04, 0.83586E+04,
     + 0.86571E+04, 0.89602E+04, 0.92680E+04, 0.95805E+04, 0.98977E+04,
     + 0.10220E+05, 0.10547E+05, 0.10878E+05, 0.11215E+05, 0.11556E+05,
     + 0.11903E+05, 0.12254E+05, 0.12611E+05, 0.12972E+05, 0.13338E+05,
     + 0.13710E+05, 0.14086E+05, 0.14468E+05, 0.14855E+05, 0.15247E+05,
     + 0.15644E+05, 0.16046E+05, 0.16453E+05, 0.16866E+05, 0.17283E+05,
     + 0.17706E+05, 0.18135E+05, 0.18568E+05, 0.19007E+05, 0.19452E+05,
     + 0.19901E+05, 0.20356E+05, 0.20817E+05, 0.21283E+05, 0.21754E+05,
     + 0.22231E+05, 0.22713E+05, 0.23201E+05, 0.23695E+05, 0.24194E+05,
     + 0.24699E+05, 0.25209E+05, 0.25725E+05, 0.26247E+05, 0.26775E+05,
     + 0.27308E+05, 0.27847E+05, 0.28393E+05, 0.28944E+05, 0.29500E+05,
     + 0.30063E+05, 0.30632E+05, 0.31207E+05, 0.31788E+05, 0.32375E+05,
     + 0.32968E+05, 0.33568E+05, 0.34173E+05, 0.34785E+05, 0.35403E+05,
     + 0.36028E+05, 0.36659E+05, 0.37296E+05, 0.37939E+05, 0.38590E+05,
     + 0.39246E+05, 0.39909E+05, 0.40579E+05, 0.41256E+05, 0.41939E+05,
     + 0.42629E+05, 0.43325E+05, 0.44029E+05, 0.44739E+05, 0.45456E+05,
     + 0.46180E+05, 0.46911E+05, 0.47649E+05, 0.48394E+05, 0.49146E+05,
     + 0.49905E+05, 0.50671E+05, 0.51445E+05, 0.52226E+05, 0.53014E+05,
     + 0.53809E+05/
  
      eps=0.01
c
      gsi = xgj(iso)
      do I=1,NT
        Q(I)=QofT(iso,I)
      enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
      SUBROUTINE AtoB(aa,bb,A,B,npt)
c***************************
c...LaGrange 3- and 4-point interpolation
c...arrays A and B are the npt data points,  given aa, a value of the 
c...A variable, the routine will find the corresponding bb value
c
c...input:  aa
c...output: bb 
      implicit DOUBLE PRECISION (a-h,o-z)
      Parameter (Nmax=119)
      dimension A(Nmax),B(Nmax)
c
C 
c
      DO 50 I=2,npt
      IF(A(I).GE.aa)THEN 
      IF(I.LT.3 .OR. I.EQ.npt) THEN
C     LaGrange three point interpolation 
      J = I
      IF(I.LT.3) J = 3
      IF(I.EQ.npT) J = npt
c.....do not devide by zero
            A0D1=A(J-2)-A(J-1)
            IF(A0D1.EQ.0.) A0D1=0.0001
            A0D2=A(J-2)-A(J)
            IF(A0D2.EQ.0.) A0D2=0.0001
            A1D1=A(J-1)-A(J-2)
            IF(A1D1.EQ.0.) A1D1=0.0001
            A1D2=A(J-1)-A(J)
            IF(A1D2.EQ.0.) A1D2=0.0001
            A2D1=A(J)-A(J-2)
            IF(A2D1.EQ.0.) A2D1=0.0001
            A2D2=A(J)-A(J-1)
            IF(A2D2.EQ.0.) A2D2=0.0001
c
      A0=(aa-A(J-1))*(aa-A(J))/(A0D1*A0D2)
      A1=(aa-A(J-2))*(aa-A(J))/(A1D1*A1D2)
      A2=(aa-A(J-2))*(aa-A(J-1))/(A2D1*A2D2)
c
      bb = A0*B(J-2) + A1*B(J-1) + A2*B(J)
c
      ELSE
C     LaGrange four point interpolation 
      J = I
c.....do not devide by zero
            A0D1=A(J-2)-A(J-1)
            IF(A0D1.EQ.0.) A0D1=0.0001
            A0D2=A(J-2)-A(J)
            IF(A0D2.EQ.0.) A0D2=0.0001
            A0D3 = (A(J-2)-A(J+1))
            IF(A0D3.EQ.0.) A0D3=0.0001
c
            A1D1=A(J-1)-A(J-2)
            IF(A1D1.EQ.0.) A1D1=0.0001
            A1D2=A(J-1)-A(J)
            IF(A1D2.EQ.0.) A1D2=0.0001
            A1D3 = A(J-1)-A(J+1)
            IF(A1D3.EQ.0.) A1D3=0.0001
c
            A2D1=A(J)-A(J-2)
            IF(A2D1.EQ.0.) A2D1=0.0001
            A2D2=A(J)-A(J-1)
            IF(A2D2.EQ.0.) A2D2=0.0001
            A2D3 = A(J)-A(J+1)
            IF(A2D3.EQ.0.) A2D3=0.0001
c
            A3D1 = A(J+1)-A(J-2)
            IF(A3D1.EQ.0.) A3D1=0.0001
            A3D2 = A(J+1)-A(J-1)
            IF(A3D2.EQ.0.) A3D2=0.0001
            A3D3 = A(J+1)-A(J)
            IF(A3D3.EQ.0.) A3D3=0.0001
c
      A0=(aa-A(J-1))*(aa-A(J))*(aa-A(J+1))
      A0=A0/(A0D1*A0D2*A0D3)
      A1=(aa-A(J-2))*(aa-A(J))*(aa-A(J+1))
      A1=A1/(A1D1*A1D2*A1D3)
      A2=(aa-A(J-2))*(aa-A(J-1))*(aa-A(J+1))
      A2=A2/(A2D1*A2D2*A2D3)
      A3=(aa-A(J-2))*(aa-A(J-1))*(aa-A(J))
      A3=A3/(A3D1*A3D2*A3D3)
c
      bb = A0*B(J-2) + A1*B(J-1) + A2*B(J) + A3*B(J+1)
      ENDIF 
c
      GO TO 100
      ENDIF 
   50 CONTINUE
  100 CONTINUE
ccc      write(2,*) 'F1, F2, F3, H1, H2, H3 =',B(J-2),B(J-1),B(J),
ccc     + A(J-2), A(J-1), A(J)
ccc      write(2,*) 'A0, A1, A2, bb =',A0,A1,A2,bb
c
      RETURN
      END 
C
C --------------------------------------------     
C*********************************************************************
      Subroutine ReadDataO2(path_to_input_files)
C*********************************************************************
C"ReadDatO2": READ data for O2
cIncluding - the CIA data (file CIAO2.dat)
C                   - thespectroscopic data for the O2 A band lines (file SDFO2.DAT)
C                   - The relaxation matrix  data (file RMFO2.DAT)
C
COutput Parameters of Routine (Arguments or Common)
C----------------------------------
Co2_nLines : Integer Array of the number of lines (Output).
C
COther important Output Quantities (through Common Statements)
C---------------------------------
C SigCIA   : WaveNumbers of CIA data file (Cm-1)
C CIA : Intensities of CIA (Cm-1Amagat-2)
C o2_Sig : WaveNumbers of the Lines (Cm-1) 
C DipoR : Rigid rotor dipole 
C Dipo0 : True Dipole transition Moments of the Lines
C E : Energies of the Lower levels of the lines (Cm-1)
C HWT0 : Air-broadened Half-Widths (at 296 K) of the 
C           Lines (Cm-1/Atm)
C BHW : Temperature Dependence Coefficients of HWT0
C o2_PopuT0 : Populations of the Lower Levels of the Lines
C                at 296 K.
C WT0 : Air-broadened Relaxation Operator Elements 
C           (at 296 K) of all Couples of Lines (Cm-1/Atm)
C BW : Temperature Dependence Coefficients of WT0
C
C Accessed Files:'SDF.dat'(Spectroscopy Data File)
C              'RMF.dat'(Relaxation Matrix File)
C 'CIA.dat'(Collision Induced Absorption File)
C
C Called By: 'ConvTP'(CONVert to Temperature and Pressure)
C---------
C     
C Double Precision Version
C     
CH. Tran last change 28 November 2005
C*********************************************************************
C
      Implicit None
      include "../gfit/ggg_int_params.f"
      Integer nLmx,nMmx,nCIAmx
      Integer iCIA,nCIA,iLine,o2_nLines,iMatrix,nMatrix
      Integer iFile,Nf,m1,lnbc
      Integer lp,idum
      character*1 R
      Double Precision SigCIA,CIA
      Double Precision o2_Sig,DipoR,Dipo0,E,HWT0,BHW,o2_PopuT0,Shift
      Double Precision WT0,BW
      Real*8 SBHW
      character path_to_input_files*(mpath+80)
C Max Number of Lines, of matrix elements 
      Parameter (nLmx=100)
      Parameter (nMmx=5000)
      Parameter (nCIAmx=100)
C Characteristic of the Band
      Common/o2_Bands/o2_nLines,nCIA
C Data of the CIA part
      Common/o2_LineSgCIA/SigCIA(nCIAmx)
      Common/SCIA/CIA(nCIAmx)
C Data of the Lines at Ref Temperature
      Common/o2_LineSg/o2_Sig(nLmx) 
      Common/DipolT0/Dipo0(nLmx) 
      Common/DipolR/DipoR(nLmx)
      Common/o2_Energy/E(nLmx) 
      Common/o2_GamT0/HWT0(nLmx) 
      Common/DTGAM/BHW(nLmx) 
      Common/DEL/SHIFT(nLmx)
      Common/o2_PopTrf/o2_PopuT0(nLmx)
      Common/NJline/m1(nLmx)
C Data of the relaxation matrix
      Common/RelxT0/WT0(nMmx)
      Common/dTRelx/BW(nMmx)

      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mauxcol  ! Avoid compiler warning (unused parameter)
      idum=mcolvav  ! Avoid compiler warning (unused parameter)
      idum=mgas     ! Avoid compiler warning (unused parameter)
      idum=mlev     ! Avoid compiler warning (unused parameter)
      idum=mrow_qc  ! Avoid compiler warning (unused parameter)
      idum=mspeci   ! Avoid compiler warning (unused parameter)
      idum=mvmode   ! Avoid compiler warning (unused parameter)
      idum=ncell    ! Avoid compiler warning (unused parameter)
      idum=nchar    ! Avoid compiler warning (unused parameter)

c File number for read
      iFile=3
C------------
C
C Open File of CIA data and Read
C
      lp=lnbc(path_to_input_files)
      iCIA=1
c     write(*,*)'path_input_to_files: ',path_to_input_files(:lp)//
c    & 'CIAO2.dat'
      Open (Unit=iFile,file=path_to_input_files(:lp)//'CIAO2.dat',
     & Status='Old')
15    Read (iFile,*,End=17)SigCIA(iCIA),CIA(iCIA)
      CIA(iCIA)=CIA(iCIA)/0.209d0  !< divide by pre-included O2 vmr
      iCIA=iCIA+1
C Check that dimensions are fine. Stop if not
      if (iCIA.GT.nCIAmx) then
         Write(*,1000)
         Stop
      Endif
1000  Format(//,1x,'************ PROBLEM !!!! ******************',
     /       /,1x,'Arrays in for CIA data storage are too small',
     /       /,1x,'raise the value of nCIAmx in ALL Parameter ',
     ,   'Statements')
C
      Goto 15
17    Close (iFile)
      nCIA=iCIA-1
C
C Open File of spectroscopic Data and Read
C        
        iLine=1
      Open(Unit=iFile,file=path_to_input_files(:lp)//'SDFO2.dat',
     &  Status='Old')
  2     Read(iFile,1001,End=100)o2_Sig(iLine),o2_PopuT0(iLine),
     , DipoR(iLine),Dipo0(iLine),E(iLine),
     , HWT0(iLine),BHW(iLine),
     , SHIFT(iLine),r,Nf
      if (r.eq.'P') m1(iLine)=-nf
      if (r.eq.'R') m1(iLine)=nf+1

c  GCT Added the following five lines Dec 11, 2011
        if(r.eq.'P')Nf=nf-1
        if(r.eq.'R')Nf=nf+1
        SBHW=0.02204D0+0.03749
     /    /(1.D0+0.05428*Nf-1.19d-3*Nf**2+2.073d-6*Nf**4)
        HWT0(iLine)=1.023*1.012*SBHW/DSQRT(1.D0+((Nf-5.D0)/55.D0)**2)

        iLine=iLine+1
C Check that dimensions are fine. Stop if not
                   If ( iLine .GT. nLmx ) Then
                   Write(*,2000)                 
                   Stop
                   End If
2000  Format(//,1x,'************ PROBLEM !!!! ******************',
     /       /,1x,'Arrays in for Line data storage are too small',
     /       /,1x,'raise the value of nLmx in ALL Parameter ',
     ,   'Statements')
C
        GoTo 2
 100  Close(iFile)
1001  Format(f12.6,2x,e9.3,2x,e9.3,1x,e9.3,2x,f10.4,2x,f5.4,2x,f4.2,
     & 2x,F8.6,4x,A1,i3,A1)
      o2_nLines=iLine-1
C        
C Open File of Relaxation Matrix Data and Read
C        
      Open(Unit=iFile,File=path_to_input_files(:lp)//'RMFO2.dat',
     & Status='Old')
      iMatrix=1
1002  Read(iFile,1003,end=1007)WT0(iMatrix),BW(iMatrix)
1003  Format(E15.7,F16.12)
c
      iMatrix=iMatrix+1
C Check that dimensions are fine. Stop if not
      If ( iMatrix .GT. nMmx ) Then
      Write(*,3000)                 
      Stop
      End If
3000  Format(//,1x,'************ PROBLEM !!!! ******************',
     /       /,1x,'Arrays in for matrix data storage are too small',
     /       /,1x,'raise the value of nMmx in ALL Parameter ',
     ,   'Statements')
C
      Goto 1002
1007  Close(iFile)
      nMatrix=iMatrix-1
      Return
      End Subroutine ReadDataO2

      SUBROUTINE RESOL(N,SYSMAT,COEF)
C
C "RESOL" : FOR RESOLUTION OF LINEAR SYSTEM
C  -----        ----          -
C
C  THIS SUBROUTINE SOLVES A LINEAR SYSTEM OF "N" UNKNONW PARAMETERS
C  BY THE GAUSS TRIANGULATION METHOD
C
C   VARIABLE
C   --------
C          N : NUMBER OF UNKNOWN QUANTITIES
C     SYSMAT : LINEAR SYSTEM MATRIX (THE N+1 COLUMN GIVES
C              THE VECTOR THAT IS ON THE RIGHT HANDSIDE OF
C              THE LINEAR SYSTEM
C       COEF : PARAMETERS THAT ARE SOLUTION OF THE LINEAR SYSTEM
C        IOK : IS EQUAL TO ZERO IF THERE IS NO PROBLEM IN SOLVING
C              THE SYSTEM. IT IS EQUAL TO 1 IF THE MATRIX IS
C              SINGULAR
C
C NO FILE ACCESSED
C NO SUBROUTINE CALLED
C CALLED BY : "ICTCVG", "ICTCVW"
C--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      DIMENSION COEF(20),SYSMAT(20,21)
      DATA EPS/1.d-30/
C
      IOK=1                                                            
      DO 300 M=1,N-1                                                    
      SUP=dABS(SYSMAT(M,M))    
      K=M
C SEARCH FOR LARGEST ELEMENT ON COLUMN
         DO 100 I=M,N
         W=dABS(SYSMAT(I,M))
         IF(W.LT.SUP) GO TO 100
         SUP=W
         K=I
100      CONTINUE
C MAKE PERMUTATION OF COLUMNS IF REQUIRED
               IF(K.NE.M)THEN
               DO 200 I=1,N+1
               BUFF=SYSMAT(K,I)
               SYSMAT(K,I)=SYSMAT(M,I)
               SYSMAT(M,I)=BUFF         
 200           CONTINUE             
               END IF
C                     QUIT COMPUTATION IF SINGULAR MATRIX
                      IF(ABS(SYSMAT(M,M)).LT.EPS) RETURN 
      DO 300 I=M+1,N                                                    
      DO 300 J=M+1,N+1                                                  
      T=SYSMAT(I,M)*SYSMAT(M,J)           
      SYSMAT(I,J)=SYSMAT(I,J)-T/SYSMAT(M,M)
300   CONTINUE                                                          
C                     QUIT COMPUTATION IF SINGULAR MATRIX
                      IF(ABS(SYSMAT(N,N)).LT.EPS) RETURN
C SOLVE TRIANGULAR SYSTEM
      COEF(N)=SYSMAT(N,N+1)/SYSMAT(N,N)     
      DO 500 I=1,N-1                                                    
      J=N-I                                                             
      SUM=0.D0                                                          
      DO 400 K=J+1,N                                                    
      SUM=SUM+SYSMAT(J,K)*COEF(K)
400   CONTINUE                                                          
      COEF(J)=(SYSMAT(J,N+1)-SUM)/SYSMAT(J,J)
500   CONTINUE                                                          
      IOK=0                                                            
      RETURN                                                            
      END  


