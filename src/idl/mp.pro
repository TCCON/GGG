      pro mp
      ncol=0
      nrow=0
      nx=0
      idir=0
      openr, unitr, 'mp.dat', /get_lun
      readf, unitr, ncol,nrow,nz,idir,frx,fry
      close, unitr
      free_lun, unitr

      for ipanel=1,ncol*nrow do begin
        mmm,ipanel,ncol,nrow,idir,frx,fry,x0,y0,x1,y1
        print, ncol,nrow,idir,x0,y0,x1,y1
      endfor
      end
