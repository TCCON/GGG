      integer function lloc(string,dl)
c  Returns last location (LLOC) of character DL in STRING
c  Useful for extracting filename from path.
c  GCT 11-Aug-2001  
c
      implicit none
      character string*(*),dl*1
      do lloc=len(string),1,-1
         if ( string(lloc:lloc).eq.dl) return
      end do
      lloc=0   ! DLM not found
      return   ! Abnormal return
      end
