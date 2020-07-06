      integer*4
     & maddln,
     & mcharhead

      character countfmt*13

      parameter(maddln=200)     ! maximum number of additional lines in a post processing output file header
      parameter(mcharhead=17200)  ! maximum number of characters in one line of a post processing output file header. Long enough for sf= line.
      parameter(countfmt='(i2,i5,i7,i4)')  ! Format of the first line in the postproc files
