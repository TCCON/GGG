pro ffile,filenm,partition,ipart,npart,path
;
;   FFILE     Version 3.1.1     GCT    24-Feb-99
;   Locates a file and returns the path to it (including the filename).
;   Systematically searches through the list of partitions, starting at
;   the partition where the previous spectrum was found.
;
;   Inputs:
;       FILENM:   C*(*)  the name of the file to locate
;
;   Outputs:
;    PARTITION:   Character array of partitions to be searched
;        IPART:   Index of current partition where last spectrum was found
;        NPART:   Number of partitions to be searched.
;         PATH:   C*(*) the full pathname to the file (including FILNAM)
;
;   Notes: 
;     1) Based on FORTRAN subroutine FINDFILE.F
;     2) PARTITION, IPART, & NPART are returned to the calling program
;        simply to save their values (no SAVE statement in IDL).
;     3) Error return:  strlen(path) = 0
;
;  Assumptions:
;     1) Relevant list of partitions is: '/ggg/mkiv/m4part.lst'
;     2) Partition names beginning with a ':' are skipped
;
;
if strlen(filenm) gt 0 then begin    ; skip empty filenames
;
; On first call, read file of data partitions to be searched.
   if npart le 0 then begin
;
;  Determine number of partitions (NPART) and define array PARTITION(NPART)
      openr,unit,'$GGGPATH/config/data_part.lst',/get_lun
      j=0
      while not EOF(unit)  do begin
         repeat  readf,unit,path  until  strmid(path,0,1) ne ':' or EOF(unit)
         j=j+1
      endwhile
      npart=j
      partition=strarr(npart)
      ipart=0
;
;  Rewind file, and read NPART partition names into array PARTITION.
      point_lun,unit,0
      for j=0,npart-1 do begin
         repeat  readf,unit,path  until  strmid(path,0,1) ne ':' or EOF(unit)
         partition(j)=path
      endfor
      close,unit
      free_lun,unit
   endif  ;  npart le 0
; 
;  Systematically search over the NPART partitions starting at IPART
   for j=0,npart-1 do begin
      path=string(partition(ipart),strtrim(filenm,2))
      print,j,path
      r=findfile(path,count=num)
      if num eq 1 then return   ; SUCCESSFUL exit
      ipart=(ipart+1) mod npart
   endfor
endif   ;   nc ge 0 
;
;   Note that after an unsucessful exit, IPART has the same value that it had
;   on entry, having been incremented NPART times, and wrapped around once.
path=''
return   ;  FAILURE EXIT  strlen(path)=0
end
