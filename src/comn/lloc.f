      integer function lloc(string,dlm)
c  Returns the Last Location Of Character DLM in STRING.
c  LLOC =  Last Location Of Character
c  GCT 15-Oct-2007   
c
      implicit none
      character string*(*),dlm*1
      do lloc=len(string),1,-1
         if( string(lloc:lloc) .eq. dlm ) return
      end do
      lloc=0   ! DLM not found anywhere
      return   ! Abnormal return
      end
