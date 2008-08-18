      subroutine clistindex(n1,list1,n2,list2,cindex)
c  Finds the locations of the strings LIST1(N1) in LIST2(N2)
c  and provides an index of the contents of LIST2 such that
c         LIST1(i) = LIST2(cindex(i))  i=1,N1
c
c  Inputs:
c     N1          I*4    Number of items in first list (LIST1)
c     LIST1(N1)   C*(*)  First Array of character strings
c     N2          I*4    Number of items in the second list (LIST2)
c     LIST2(N2)   C*(*)  Second array of character strings
c
c  Outputs:
c     CINDEX(N1)  I*4    Indexing array.
c
c  Notes:
c  1) Output array (CINDEX) has a dimension of N1 and a maximum value of N2.
c  2) If CINDEX(i) > 0  LIST1(i) was found in LIST2 at location CINDEX(i)
c  3) if LIST1(i) was not found anywhere in LIST2, then CINDEX(i) is unchanged
c  4) If LIST2 contains duplicate entries, only the first one will be indexed.
c  5) Array CINDEX should normally be set to zero before calling CLISTINDEX.
c  6) CLISTINDEX is case sensitive. To make it case insensitive you should
c     call LOWERCASE(LIST1) and/or LOWERCASE(LIST1) prior to calling CLISTINDEX.
c  7) N1 can be >, =, or < N2.
c  8) If N1=1, then CLISTINDEX simply searches for the single string LIST1
c     in the string array LIST2. Its location is returned in CINDEX, which
c     is a simple I*4 integer, rather than an array.
c
c Example 1:
c LIST1 = h2o h2o h2o co2 o3  n2o co  ch4 
c LIST2 = co2 h2o ch4 hcl h2o 
c CINDEX=  2   2   2   1   0   0   0   3

c Example 2:
c LIST1 = co2 h2o ch4 hcl h2o 
c LIST2 = h2o h2o h2o co2 o3  n2o co  ch4 
c CINDEX = 4   1   8   0   1 
c
      implicit none
      integer*4 n1,i1,n2,i2,cindex(n1)
      character list1(n1)*(*),list2(n2)*(*)
c
      do i1=1,n1
         if(cindex(i1).eq.0) then
            do i2=n2,1,-1
               if(list1(i1).eq.list2(i2)) cindex(i1)=i2
            end do
         endif
      end do
      return
      end
