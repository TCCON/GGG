      function ciso2kiso(kgas,ciso)
      integer*4 kgas,ciso2kiso
      character*1 ciso

      ciso2kiso=ichar(ciso)-48
      if(kgas.eq.2) then
         if(ciso2kiso.ge.17) then
            ciso2kiso=ciso2kiso-6   ! Isotope A -> 11; B -> 12
         elseif(ciso2kiso.eq.0) then
            ciso2kiso=10       ! Isotope  0 -> 10
         elseif(ciso2kiso.gt.9) then
            write(*,*) 'kgas,ciso = ',kgas,ciso
            stop 'ciso2kiso: Unknown CO2 isotopolog'
         endif
      endif
      return
      end

