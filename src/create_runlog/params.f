
      integer*4
     & bytepw,                 !Number of bytes per data word (=2 for AT,M4; =4 for OP,GR)
     & ifirst,                 !The index of the first spectral point on disk
     & ilast,                  !The index of the last spectral point on disk
     & possp,                  !Number of bytes before IFIRST'th point (i.e. header length)
     & object,                 !Heavenly object (Moon=1, Sun=2)
     & ios,                    !Status flag (0=success, 1=EOF)
     & iy,im,id,jd

