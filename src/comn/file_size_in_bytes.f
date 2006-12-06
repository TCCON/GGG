       function file_size_in_bytes(lunr,file_path)
c  Returns the size of any file (in bytes)
c
c  Inputs:
c      LUNR               I*4  Logical Unit Number
c      FILE_PATH          C*(*) Full path name
c
c  Output:
c      file_size_in_bytes I*4  File size (bytes)
c
c  Performs direct-access reads to detect the EOF location.
c  Keeps track of which reads are within the file and which
c  fall off the end.
c  For a large file such as the atm.101 linelist (117 Mbytes)
c  It will take 27 reads to fall off the EOF and then another
c  27 reads to determine the exact EOF position. This is much
c  faster than reading the entire file from beginning to end.
c
       integer*4 file_size_in_bytes, lunr, new, nhi, nlo, istat
       character file_path*(*),string*1
       open(lunr,file=file_path,access='direct',
     $ form='unformatted',status='old',recl=1)
c
c Initialize
       nlo=0
       nhi=0
       istat=0
c
c  Keep doubling NHI until it falls off the end of file (EOF)
       do while (istat.eq.0)
          nlo=nhi
          nhi=2*nhi+1
          read(lunr,rec=nhi,iostat=istat) string
       end do
c  So now we know that the EOF lies between NLO and NHI
c
c  Find the exact EOF location by iterative bisection.
       do while (nhi-nlo.gt.1)
          new=(nhi+nlo)/2
          read(lunr,rec=new,iostat=istat) string
          if(istat.eq.0) then
            nlo=new
          else
            nhi=new
          endif
        end do
c
       close(lunr)
       file_size_in_bytes=nlo
       return
       end
