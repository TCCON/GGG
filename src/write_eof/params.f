
      integer*4 lun_qc, mcol, mrow, mluns
      integer*4 lun_mul, lun_rlg

      parameter (lun_mul=21)       ! multiggg.sh 
      parameter (lun_rlg=22)       ! runlog
      parameter (lun_qc=26)
      parameter (mcol=150)
      parameter (mrow=9999)
      parameter (mluns=100)

      character
     & gfit_version*80,
     & gsetup_version*80,
     & atmsum*32,
     & gctsum*32,
     & fciasum*32,
     & sciasum*32,
     & solarsum*32,
     & csformat*105,
     & header*22000
