      integer function posnall (unit,fpos,nrec)
c
c  Version:   2.1.0    16-5-95     GCT
c
c  Description:
c       Positions in an bounded (i.e. known length) HITRAN-format ascii
c       linelist to the last line having a frequency <= FPOS.
c
c  On Input:
c       unit  = logical unit number of already-opened linelist
c       fpos  = frequency (cm-1) that you want to position to
c       nrec  = number of records (lines) in the file
c
c  On Output:
c       posnall  = record # of last line whose freq < fpos
c                = 0 means that all lines exceeded fpos
c                = NREC means that no line exceeded fpos
c
c  Normal usage:
c     fsib=file_size_in_bytes(lun_ll,linfile)
c     nlines=fsib/reclen
c     open(lun_ll,file=linfile,access='direct',
c    & form='formatted',status='old',recl=reclen)
c     k1=posnall(lun_ll,nu1,nlines) ! index of the last line with v < NU1
c     k2=posnall(lun_ll,nu2,nlines) ! index of the last line with v < NU2
c     do kk=k1+1,k2
c        read(lun_ll,llformat,rec=kk) igas, iso, freq, strength, eprime, ....
c     end do

      implicit none
      integer unit,nrec,new,nlo,nhi
      real*8 freq,fpos
C======================================================================
c  Position to frequency FPOS by iterative bisection (divide and conquer),
c  using a variant of the LOCATE subroutine described in Numerical Recipes.
c  The total number of iterations required is roughly LOG2(NLINE),
c  which for the HITRAN list (10^6 lines) is 20.
      nlo=0                         ! don't worry, NEW can never be < 0
      nhi=nrec+1                    ! NHI can never exceed nrec
      do while(nhi-nlo.gt.1)
         new=(nlo+nhi)/2            ! Bisect remaining range of positions
         read(unit,'(3x,f12.6)',rec=new) freq
         if(freq.lt.fpos) then      ! FPOS lies between records NLO and NHI
            nlo=new
         else
            nhi=new
         endif
      end do
      posnall=nlo
      return
      end
