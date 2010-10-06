      SUBROUTINE GETLUN ( UNIT )
 
      integer unit

      CALL FNDLUN ( UNIT )
      
      if(unit.eq.0) return
      
      call reslun(unit)

      end subroutine getlun

      subroutine freelun(unit)
      
      integer unit
      integer iostat
      logical opened

      if (unit.gt.0) then
         INQUIRE ( UNIT=unit, OPENED=OPENED, IOSTAT=IOSTAT )
         if(opened) then
            close(unit)
         else
            write(*,*) "Tried to close unoped file"
         endif
         
         call frelun(unit)
      else
         write(*,*) "Failed closing: ", unit
      endif

      end subroutine freelun



