      subroutine read_mav_head(lunm,nlev,nspeci,nlhead)
c  Reads the contents of the .mav file MAVNAME header into the appropriate arrays.
c
c  Inputs:
c       LUNM           Logical unit number
c  Outputs:
c       NLEV           Number of levels of MAVFILE to be read
c       NLHEAD         Number of header lines
c       NSPECI         Actual number of different gases in MAVFILE
c
      implicit none
      integer*4 lunm,nlev,nspeci,nlhead,ncol

      read(lunm,*)nlhead,ncol,nlev
      nspeci = ncol-4
      return
      end
