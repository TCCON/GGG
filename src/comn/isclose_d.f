      function isclose_d(a,b)
c Test if two double-precision values (a & b) are 
c close enough in value to be considered equal. 
c Specifically, they are considered close enough if
c
c   |a - b| < max(|a|, |b|) * 2^-51
c
c Inputs
c
c   A, B : REAL*8 - values to compare
c
c Returns
c
c   .TRUE. is A and B are close in value, .FALSE. otherwise

      real*8 a, b, eps
      logical isclose_d

c In double-precision (real*8) floats, the non-exponent value of
c the number uses 52 bits. That means the spacing between two real*4s
c with the same exponent is 2**-52. However in the worst case, A
c might round down and B round up (or vice verse) so we allow 2x this
c spacing, or 2**-51. 
c
c To account for the magnitude, we scale by the larger of the two
c values.
      eps = max(abs(a), abs(b)) * 2.d0**(-51.0d0)
c Must be less than or equal to to work if both numbers are 0
      isclose_d = abs(a - b) .le. eps 

      return
      end
