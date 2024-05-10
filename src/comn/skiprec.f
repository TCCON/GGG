      subroutine skiprec(lun,nrec)
c  Skips NREC records of the already open logical unit number LUN.
c  Avoids having multiple read(lun,*) or backspace(lun) statements in
c  the main program or setting up do-loops containing read/backspaces.
      implicit none
      integer*4 lun,irec,nrec
      if (nrec.gt.0) then
         do irec=1,nrec
            read(lun,*)
         end do
      elseif(nrec.lt.0) then
         do irec=nrec,-1
            backspace(lun)
         end do
      endif
      return
      end
