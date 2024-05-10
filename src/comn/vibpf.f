      function vibpf(temp,vibfrq,dgen,nvmode)
c Calculates the vibrational partition function:
c the fraction of molecules that are in the vibrational ground state.
      integer*4 i,nvmode,dgen(nvmode)
      real*4 temp,vibfrq(nvmode),vibpf
      vibpf=1.0
      do i=1,nvmode
         vibpf=vibpf*(1.-exp(-(1.43881*vibfrq(i)/temp)))**dgen(i)
      end do
      return
      end
