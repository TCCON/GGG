      subroutine default_run_info(mif,mns,infomat)
c
c  Input:
c    mif       I*4    Maximum number of info channels
c    mns       I*4    Maximum number of interferograms per scan set (max NSS)
c
c  Output:
c    infomat(mif,mns)R*8 Information produced by real-time algorithm
c
      implicit none

      integer*4
     & mif,        ! Subroutine input argument (see above)
     & mns,        ! Subroutine input argument (see above)
     & scanind     ! Scan index

      real*8
     & infomat(mif,mns)! Subroutine output argument (see above)
c
c  Fill output vectors with default or "not found" values.
      do scanind=1,mns
         infomat(1,scanind)=100.0d0       ! HUM  IFHum
         infomat(2,scanind)=27.0d0        ! TLP  IFSSrcT
         infomat(3,scanind)=0.71972d0     ! PIM  IFS_P
         infomat(4,scanind)=28.1d0        ! TSC  ScBlkl_T
         infomat(5,scanind)=1500.0d0      ! PGR  InGaAs_R
         infomat(6,scanind)=2800.0d0      ! PGR  Si_R
         infomat(7,scanind)=45.9448d0     ! LAT  Latitude
         infomat(8,scanind)=-90.2732d0    ! LON  Longitude
         infomat(9,scanind)=442.0d0       ! ALT  Altitude
         infomat(10,scanind)=2.7d0        ! WSA  Zeno_WindSpeed_avg
         infomat(11,scanind)=2.2d0        ! WSS  Zeno_WindSpeed_std
         infomat(12,scanind)=7.9d0        ! WSM  Zeno_WindSpeed_max
         infomat(13,scanind)=106.0d0      ! WDA  Zeno_WindDir_avg
         infomat(14,scanind)=79.6d0       ! WDS  Zeno_WindDir_std
         infomat(15,scanind)=24.3d0       ! TOU  Zeno_Temp_avg
         infomat(16,scanind)=48.5d0       ! HOU  Zeno_RH_avg
         infomat(17,scanind)=636.7d0      ! ZSA  Zeno_SolarRadiance_avg
         infomat(18,scanind)=0.8d0        ! ZSS  Zeno_SolarRadiance_std
         infomat(19,scanind)=965.1d0      ! POU  Zeno_Press_avg
         infomat(20,scanind)=0.0d0        ! ZRM  Zeno_Rain_max
         infomat(21,scanind)=9.0d0        ! ZLM  Zeno_Lightning_max
         infomat(22,scanind)=15.1d0       ! ZVA  Zeno_VBatt_avg
         infomat(23,scanind)=206.9d0      ! DAA  Dome_azi_avg
         infomat(24,scanind)=0.0d0        ! DSM  Dome_Status_max
         infomat(25,scanind)=205.7d0      ! SAA  ST_tpg_azi_avg
         infomat(26,scanind)=42.1d0       ! SEA  ST_tpg_ele_avg
         infomat(27,scanind)=0.0d0        ! STM  ST_TPS_max
         infomat(28,scanind)=-1.0d0       ! SIA  ST_t_int_avg  VALUE IMPORTANT
         infomat(29,scanind)=-1.0d0       ! SIS  ST_t_int_std  VALUE IMPORTANT
         infomat(30,scanind)=-0.5d0       ! SOA  ST_off_azi_avg
         infomat(31,scanind)=-0.1d0       ! SOE  ST_off_ele_avg
         infomat(32,scanind)=3.0d0        ! SDA  ST_Tdrift_avg
         infomat(33,scanind)=-9999999.0d0 ! IDA  IFSDT_avg IFS125_TIME - NTP_TIME
         infomat(34,scanind)=0.0d0        ! SFM  ScanType
         infomat(35,scanind)=-1.0d0       ! ISS  ScanStatus    VALUE IMPORTANT
      enddo
      return
      end
