c
c  This is the list of indices used with the header I*4 vector
c
      integer*4
     & i_ssm,      ! Index of Sample Spacing Multiplier
     & i_ssp,      ! Index of Sample Spacing divisor
     & i_srccode,  ! Index of numerical code of light source selection
     & i_sfmcode,  ! Index of numerical code of scan type (Sample ForM)
     & i_dtccode,  ! Index of numerical code of detector
     & i_pgn,      ! Index of pre-amplifier gain master
     & i_sgna,     ! Index of main amplifier binary gain slave
     & i_sgnb,     ! Index of main amplifier binary gain master
     & i_hpf,      ! Index of setting of the electronics high-pass filter
     & i_gfw,      ! Index of count of Good ForWard scans
     & i_gbw,      ! Index of count of Good BackWard scans
     & i_bfw,      ! Index of count of Bad ForWard scans
     & i_bbw,      ! Index of count of Bad BackWard scans
     & i_rsn,      ! Index of Running Sample Number
     & i_pka,      ! Index of PeaK Amplitude of the forward run master chan.
     & i_pkl,      ! Index of PeaK Location of the forward run master chan.
     & i_pra,      ! Index of Peak Amplitude of the Reverse run master chan.
     & i_prl,      ! Index of Peak Location of the Reverse run master chan.
     & i_p2a,      ! Index of Peak Amplitude of the forward run slave chan.
     & i_p2l,      ! Index of Peak Location of the forward run slave chan.
     & i_p2r,      ! Index of Peak amplitude of the Reverse run slave chan.
     & i_p2k,      ! Index of PeaK location of the reverse run slave chan.
     & i_zff,      ! Index of Zero Filling Factor
     & i_inscode,  ! Index of numerical code of INStrument identifier
     & i_bmscode,  ! Index of numerical code of BeaM Splitter
     & i_aqmcode   ! Index of numerical code of AQuisition Mode
      parameter (i_ssm    = 1)
      parameter (i_ssp    = 2)
      parameter (i_srccode= 3)
      parameter (i_sfmcode= 4)
      parameter (i_dtccode= 5)
      parameter (i_pgn    = 6)
      parameter (i_sgna   = 7)
      parameter (i_sgnb   = 8)
      parameter (i_hpf    = 9)
      parameter (i_gfw    =10)
      parameter (i_gbw    =11)
      parameter (i_bfw    =12)
      parameter (i_bbw    =13)
      parameter (i_rsn    =14)
      parameter (i_pka    =15)
      parameter (i_pkl    =16)
      parameter (i_pra    =17)
      parameter (i_prl    =18)
      parameter (i_p2a    =19)
      parameter (i_p2l    =20)
      parameter (i_p2r    =21)
      parameter (i_p2k    =22)
      parameter (i_zff    =23)
      parameter (i_inscode=24)
      parameter (i_bmscode=25)
      parameter (i_aqmcode=26)
c
c  This is the list of indices used with the header R*8 vector
c
      integer*4
     & i_laserate, ! Index of laser fringe rate in Hz (1000*vel from OPUS)
     & i_lwn,      ! Index of Laser WaveNumber in cm^-1
     & i_foc,      ! Index of collimator focal length in mm
     & i_aptval,   ! Index of linear diameter of field stop in mm
     & i_lpf,      ! Index of corner frequ. of the electronics low-pass filter
     & i_dur,      ! Index of measurement duration
     & i_mvd       ! Index of Maximum Velocity Deviation
      parameter (i_laserate= 1)
      parameter (i_lwn     = 2)
      parameter (i_foc     = 3)
      parameter (i_aptval  = 4)
      parameter (i_lpf     = 5)
      parameter (i_dur     = 6)
      parameter (i_mvd     = 7)
c
c  We define non-Bruker header items, jokingly called TLA for three-letter
c  acronyms.  Our file creation subroutine can be limited to produce none,
c  a few, or all of our specific items.
c
c GCT      integer*4
c GCT     & tla_ext_none,! Use only the Bruker-defined header items
c GCT     & tla_ext_min, ! Add those that apply to all TCCON sites
c GCT     & tla_ext_pgr, ! Also add the Pre-amp Gain Resistor value
c GCT     & tla_ext_full ! Complete set from automated system (Zeno, ST, dome,...)
c GCT      parameter (tla_ext_none = 0)
c GCT      parameter (tla_ext_min  = 1)
c GCT      parameter (tla_ext_pgr  = 2)
c GCT      parameter (tla_ext_full = 3)
c
c  Header items that are specified with strings get converted to integer
c  code values.  These are easier to handle inside the main sections of
c  the software (e.g. for comparisons).  Multiple input strings can easily
c  be mapped to the same code value.  When creating OPUS files, these
c  code values get converted back to strings which can be different from
c  their input counterparts.
c
      integer*4
     & src_sun,     ! Code for source set to solar (aka emission) port
     & src_nir,     ! Code for source set to Near Infra-Red (tungsten lamp)
     & src_mir,     ! Code for source set to Mir Infra-Red (globar)
     & src_off,     ! Code for source turned off
     & sfm_solar,   ! Code for sample form of solar scan
     & sfm_cell,    ! Code for sample form of cell run
     & sfm_idle,    ! Code for sample form of idle scan
     & sfm_script,  ! Code for ancient sample form of 'scripting'
     & sfm_aeros,   ! Code for sample form of aerosol scan
c     & dtc_ingaas,  ! Code for detector of InGaAs
c     & dtc_si,      ! Code for detector of Si
c     & dtc_insb,    ! Code for detector of InSb
c     & dtc_hgcdte,  ! Code for detector of HgCdTe
c     & dtc_bolom,   ! Code for detector of bolometer
c     & ins_120,     ! Code for instrument of IFS120HR
c     & ins_120m,    ! Code for instrument of IFS120M
c     & ins_125,     ! Code for instrument of IFS125HR
c     & ins_66s,     ! Code for instrument of IFS66/S
     & bms_caf2,    ! Code for beam splitter of CaF2
     & bms_kbr,     ! Code for beam splitter of KBr
     & bms_quartz,  ! Code for beam splitter of Quartz
     & bms_sica,    ! Code for beam splitter of SiCa
     & aqm_sd,      ! Code for aquisition mode of forward + reverse
     & aqm_sf,      ! Code for aquisition mode of forward only, fast return
     & aqm_sn       ! Code for aquisition mode of forward only
      parameter (src_sun    = 1)
      parameter (src_nir    = 2)
      parameter (src_off    = 3)
      parameter (src_mir    = 4)
      parameter (sfm_solar  = 1)
      parameter (sfm_cell   = 2)
      parameter (sfm_idle   = 3)
      parameter (sfm_script = 4)
      parameter (sfm_aeros  = 5)
c      parameter (dtc_ingaas = 1)
c      parameter (dtc_si     = 2)
c      parameter (dtc_insb   = 3)
c      parameter (dtc_hgcdte = 4)
c      parameter (dtc_bolom  = 5)
c      parameter (ins_120    = 1)
c      parameter (ins_120m   = 2)
c      parameter (ins_125    = 3)
c      parameter (ins_66s    = 4)
      parameter (bms_caf2   = 1)
      parameter (bms_kbr    = 2)
      parameter (bms_quartz = 3)
      parameter (bms_sica   = 4)
      parameter (aqm_sd     = 1)
      parameter (aqm_sf     = 2)
      parameter (aqm_sn     = 3)
c
