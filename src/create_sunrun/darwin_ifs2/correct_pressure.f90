    subroutine correct_pressure(specname,gggdir,iy,im,id,hh,mm,pout)

!	Correct Zeno pressure, version Oct 2015  DG/VV
!	From 2007, add +0.8 hPa to Zeno    
!	From Nov 2014 - June 2015 Zeno power supply became noisy, pressures noisy and low offset
!	Rob Chapman provided BoM 30 min station and mean sea level pressures, file DarwinPressure_1410_1506.txt
!	For the dodgy period, use the BoM pressures adjusted to TCCON altitude with correction to be consistent with previous zeno offset.

!	Inputs/outputs
!		specname 			- required only for debug output
!		gggdir 				- GGGPATH required to locate BoM lookup file DarwinPressure_1410_1506.txt
!		iy, im, id, hh, mm 	- Date and time of spectrum (UT)
!		pout 				- on input, raw zeno pressure from the spectrum file header
!		pout 				- on output, corrected pressure for sunrun

    implicit none
    integer(4)::    ierr
    integer(4)::    iy, im, id, hh, mm
    integer(4)::    iy2, im2, id2, hh2, mm2
    integer(4)::    julian					!Function version of Julian.f, see end of this .f90 file
    real(8)::       pout, pbom, LastPbom, pout_cor, pzeno, pbom_cor
    real(8)::       SpecTimeUT, BomTimeUT, LastBomT
	character*(*)	specname, gggdir
    character*20 	blank

	blank = "BoM                 "
    pout_cor = 0.8 ! mbar based on nmd03 email: 2007-03-05
    pbom_cor = 1.0 ! correction to Pbom to fit historical Zeno corrected pressure
	
!   Begin
	pzeno = pout
	pout = pout + pout_cor		!Default correction
	
!	Correction Oct 2014 - Jun 2015 for faulty zeno battery.  Use BoM pressures
    if(julian(iy,im,id).ge.julian(2014,09,30).and.julian(iy,im,id).le.julian(2015,06,30))then
		SpecTimeUT = dble(julian(iy,im,id)) + ((dble(hh))*3600. + dble(mm)*60.)/86400.
		print *, trim(specname), iy, im, id, hh, mm, pout
		open (999, file=trim(gggdir)//"src/create_sunrun/darwin_ifs2/DarwinPressure_1410_1506.txt")
		read(999,*)
		do
			read (unit=999, fmt='(10x,I2,x,I2,x,I4,x,i2,x,i2,x,F6.1)', iostat=ierr) id2, im2, iy2, hh2, mm2, pbom
			if(ierr<0)then
				stop "Reached end of BoM pressure file"
			endif
			pbom = pbom * exp(-32./8500.)	!Correct msl pressure to 32 m
			BomTimeUT = dble(julian(iy2,im2,id2)) + ((dble(hh2)-9.5)*3600. + dble(mm2)*60.)/86400.
			if (BomTimeUT < SpecTimeUT) then
				LastPbom = pbom
				LastBomT  = BomTimeUT
			else
				pout = ((SpecTimeUT-LastBomT)*pbom + (BomTimeUT-SpecTimeUT)*LastPbom)/(BomTimeUT - LastBomT) + pbom_cor
				print *, blank,iy2, im2, id2, hh2, mm2, pout, pzeno+pout_cor-pout
				print *
				exit
			endif         
		enddo
	endif
	close(999)
	return
	end        

	function julian(y,m,d)
!	Function converts Gregorian calendar date
!	to Julian Day number at Greenwich Mean Noon
    implicit none
    INTEGER*4 Y,M,D, julian
	
!******************************************************
!        Input:
!     Y            Integral calendar year
!     M            Integral calendar month
!     D            Integral calendar day
!
!     Output
!     Julian       Integral Julian Day number

!******************************************************
     julian=367*y-7*(y+(m+9)/12)/4-3*((y+(m-9)/7)/100+1)/4+275*m/9+d+1721029
	 
     return
     end
