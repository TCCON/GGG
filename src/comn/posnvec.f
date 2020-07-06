      integer function posnvec(nrec,freq,fpos)
c
c  Finds the index of the last line having freq(i) <= FPOS
c  in an ascending, ordered list.
c
c  Inputs:
c       nrec  = number of lines in the freq-vector
c       freq(nrec) = frequencies of spectral lines
c       fpos  = frequency (cm-1) that you want to position to
c
c  Output:
c       posnvec  = record # of last line freq(posnvec) < fpos
c                = 0 means that all lines exceeded fpos
c                = NREC means that no line exceeded fpos
c
c  Positions to frequency FPOS by iterative bisection (divide and conquer),
c  using a variant of the LOCATE subroutine described in Numerical Recipes.
c  The total number of iterations required is roughly LOG2(NLINE),
c  which for the HITRAN list (8x10^6 lines) is 23.

      implicit none
      integer nrec,new,nlo,nhi
      real*8 freq(nrec),fpos

      nlo=0                         ! don't worry, NEW can never be < 0
      nhi=nrec+1                    ! NHI can never exceed nrec
      do while(nhi-nlo.gt.1)
         new=(nlo+nhi)/2            ! Bisect remaining range of positions
         if(freq(new).lt.fpos) then      ! FPOS lies between records NLO and NHI
            nlo=new
         else
            nhi=new
         endif
      end do
      posnvec=nlo
      return
      end
