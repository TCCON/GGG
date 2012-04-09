C
C Parameters to be included in the program
C
	Integer*4 nSigmx,iFile
	Integer*4 nBmx,nLmx,nIsotp
      	Integer*4 Jmax,Nlifmax,Nwmax

C Max number of spectral points
      Parameter (nSigmx=500000)
C Unit Number for File access
      Parameter (iFile=3)
C Max Number of Bands
      Parameter (nBmx=600)	
C Max Number of Lines per Band
      Parameter (nLmx=400)
C Number of CO2 Isotopomers 
      Parameter (nIsotp=8)
c Max number of l values
      parameter(Nlifmax=7)
c Max number of j values
      parameter(Jmax=125)
c Max Number of W coupling
      parameter(Nwmax=100000)

	real*8 Ct,T0,CtGamD,aMolAtm,Pi
      	real*8 aMass(nIsotp)
C Constants      
      Data Ct/1.4387686d0/
      Data T0/296.d0/
      Data CtGamD/1.1325d-08/
      Data aMolAtm/7.33889d+21/
      Data Pi/3.141592654d0/
      Data aMass/44.D-3,45.d-3,46.d-3,45.d-3
     ,     ,47.d-3,46.d-3,48.d-3,47.d-3/

C      
C Polynomial Coefficients for the Various Isotopomers
	integer*4 Jcoef
      	real*8 co2_Qcoef(nIsotp,4)
C Iso=1, i.e. O(16)-C(12)-O(16)
      DATA (co2_Qcoef(1,Jcoef),Jcoef=1,4)/-.13617D+01, .94899D+00,
     +               -.69259D-03, .25974D-05/
C Iso=2, i.e. O(16)-C(13)-O(16)
      DATA (co2_Qcoef(2,Jcoef),Jcoef=1,4)/-.20631D+01, .18873D+01,
     +               -.13669D-02, .54032D-05/
C Iso=3, i.e. O(16)-C(12)-O(18)
      DATA (co2_Qcoef(3,Jcoef),Jcoef=1,4)/-.29175D+01, .20114D+01,
     +               -.14786D-02, .55941D-05/
C Iso=4, i.e. O(16)-C(12)-O(17)
      DATA (co2_Qcoef(4,Jcoef),Jcoef=1,4)/-.16558D+02, .11733D+02,
     +               -.85844D-02, .32379D-04/
C Iso=5, i.e. O(16)-C(13)-O(18)
      DATA (co2_Qcoef(5,Jcoef),Jcoef=1,4)/-.44685D+01, .40330D+01,
     +               -.29590D-02, .11770D-04/
C Iso=6, i.e. O(16)-C(13)-O(17)
      DATA (co2_Qcoef(6,Jcoef),Jcoef=1,4)/-.26263D+02, .23350D+02,
     +               -.17032D-01, .67532D-04/
C Iso=7, i.e. O(18)-C(12)-O(18)
      DATA (co2_Qcoef(7,Jcoef),Jcoef=1,4)/-.14811D+01, .10667D+01,
     +               -.78758D-03, .30133D-05/
C Iso=8, i.e. O(17)-C(12)-O(18)
      DATA (co2_Qcoef(8,Jcoef),Jcoef=1,4)/-.17600D+02, .12445D+02,
     +               -.91837D-02, .34915D-04/
c...   CO2  --   626
c...   CO2  --   636
c...   CO2  --   628
c...   CO2  --   627
c...   CO2  --   638
c...   CO2  --   637
c...   CO2  --   828
c...   CO2  --   728
