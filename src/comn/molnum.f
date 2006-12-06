      integer function molnum(name)
c  Finds the molecule number (MOLNUM) given a left-justified string (NAME)
c  This is accomplished by searching through array 'LIST' which
c  contains a list of the molecule numbers and permissible names.
c  A case-insensitive test is performed by first converting NAME to upper case.
c  MOLNUM=0 is returned if the molecule in NAME could not be identified.
c
      integer i,jn,igas,ngas
      parameter (ngas=79)
      character list(ngas)*16,name*(*),ucname*16
c  USE UPPER CASE ONLY !!!!!!
      data list/
     &'01 H2O',                                                     ! 1
     &'02 CO2','02 WCO2','02 ICO2','02 SCO2','02 VCO2','02 TCO2',   ! 7
     &'03 O3','04 N2O','05 CO','06 CH4','07 O2','08 NO','09 SO2',   ! 14
     &'10 NO2','11 NH3','12 HNO3','13 OH','14 HF','15 HCL',         ! 20
     &'16 HBR','17 HI','18 CLO','19 OCS','20 H2CO',                 ! 25
     &'21 HOCL','22 HO2','23 H2O2','24 HONO','25 HO2NO2','25 HNO4', ! 31
     &'26 N2O5','27 CLNO3','28 HCN','29 CH3F','30 CH3CL','31 CF4',  ! 37
     &'32 CF2CL2','32 CCL2F2','32 F12',                             ! 40
     &'33 CFCL3','33 CCL3F','33 F11',                               ! 43
     &'34 CH3CCL3','35 CCL4',                                       ! 45
     &'36 CF2O','36 COF2',                                          ! 47
     &'37 CFCLO','37 COFCL','38 C2H6','39 C2H4','40 C2H2','41 N2',  ! 53
     &'42 CHF2CL','42 F22',                                         ! 55
     &'43 COCL2','44 CH3BR','45 CH3I','46 HCOOH','47 H2S',          ! 60
     &'48 F21','48 CHFCL2',                                         ! 62
     &'49 HDO','50 SF6',                                            ! 64
     &'51 F113','51 CF3CCL3',                                       ! 66
     &'52 CLCN','53 F142B',                                         ! 68
     &'54 O3-668','55 O3-668','56 O3-668',                          ! 71
     &'57 O3-1','58 O3-2','59 O3-3',                                ! 74
     &'60 HCL-1','61 HCL-2','62 HCL-3','63 HCL-4','64 NO3'/         ! 79
c
c  Convert NAME to upper case (UCNAME).  [GAS is already upper case]
      ucname=name
      do i=1,len(name)
        jn=ichar(name(i:i))
        if(jn.ge.92) ucname(i:i)=char(jn-32)  ! convert to upper case
      end do 
c
c  Search through list for matching name
      do igas=1,ngas
      if(ucname.eq.list(igas)(4:)) then    !  success
      molnum=10*(ichar(list(igas)(1:1))-48)+ichar(list(igas)(2:2))-48
      return
      endif
      end do
      molnum=0    ! It did not find a match
      return
      end
