      subroutine write_binary_file
     & (lunw,outname,header,headlen,bytepw,r4buf,nmp)
c  writes a binary file consisting of a header of length headlen
c  plus an NMP-point data vector. Depending on the value of bytepw
c  the data are written in a R*4 or I*2 format. In the latter case,
c  the data values are multiplied by 20000 before conversion to I*2,
c  under the assumption that the largest data value will be around 1.
c
      integer*4 nmp,lunw,i,headlen,bytepw,iabpw
      real*4 r4buf(nmp)
      character header*(*),outname*(*)

      iabpw=iabs(bytepw)
      open(lunw,file=outname,access='direct', status='unknown',
     & form='unformatted',recl=headlen+iabpw*nmp)
      if(iabpw.eq.2) then
         write(lunw,rec=1) header(1:headlen),
     &   (int2(r4buf(i)),i=1,nmp)
      elseif(iabpw.eq.4) then
         write(lunw,rec=1) header(1:headlen),(r4buf(i),i=1,nmp)
      else
         write(*,*)'BYTEPW=',bytepw
         stop 'Error in write_binary_file: Unknown format'
      endif
      close(lunw)
      return
      end

