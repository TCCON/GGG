      real*4 function twodint(vali,valj,tabi,tabj,tab,mi,ni,nj)
c  Performs bi-linear (2-D) interpolation into a table (TAB) of data values.
c  First searches vectors (TABI, TABJ) of I and J ordinates to bracket requested
c  point. Ordinate vectors must be sorted (ascending or descending) but need not
c  be equally spaced. Then performs bi-linear interpolation in table (TAB).
c  Returns edge value if requested point lies outside table.
c
c  Input Paramaters:
c        VALI, VALJ           co-ordinates of the requested point
c        TABI(NI), TABJ(NJ)   vector of ordinates of the tabulated points 
c        TAB(MI,NJ)           matrix of tabulated values
c        MI, NI               number of rows of matrix (MI=declared; NI=actual)
c        NJ                   number of actual columns of matrix 
c
      implicit none
      integer mi,ni,i,nj,j
      real*4 vali,valj,fi,fj,omfi,omfj,tabj(nj),tabi(mi),tab(mi,nj)

      fi=0.0 ! prevent compiler warning (may be uninitialized)
      fj=0.0 ! prevent compiler warning (may be uninitialized)
c
c  Search for TABI's which bracket VALI and find its fractional position (FI)
      do i=1,ni-1
        fi=(vali-tabi(i))/(tabi(i+1)-tabi(i))
        if(fi.ge.0.0 .and. fi.le. 1.0) go to 33   !  Successful bracketing
      end do
      if(fi.lt.0) then      ! Outside table: set to uppermost row
        i=1
        fi=0.0
      elseif(fi.gt.1) then  ! Outside table: set to lowest row
        i=ni-1
        fi=1.0
      else
        stop ' TWODINT: this should never happen. Notify GCT'
      endif
33    continue
c
c  Search for TABJ's which bracket VALJ and find its fractional position (FJ)
      do j=1,nj-1
        fj=(valj-tabj(j))/(tabj(j+1)-tabj(j))
        if(fj.ge.0.0 .and. fj.le. 1.0) go to 22   !  Successful bracketing
      end do
      if(fj.lt.0) then       ! Outside table: set to leftmost column
        j=1
        fj=0.0
      elseif(fj.gt.1) then   ! Outside table: set to rightmost column
        j=nj-1
        fj=1.0
      else
        stop ' TWODINT: this should never happen'
      endif
22    continue
c
c  Perform bi-linear interpolation
      omfi=1-fi
      omfj=1-fj
      twodint = tab(i,j)*omfi*omfj + tab(i,j+1)*omfi*fj
     &        + tab(i+1,j)*fi*omfj + tab(i+1,j+1)*fi*fj 
      return
      end
