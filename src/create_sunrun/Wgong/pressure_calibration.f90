    subroutine pressure_calibration (Pout, specname, pout_corr)

!   This subroutine corrects the raw measured pressure following
!   post-measurement calibration of the pressure sensor.
!   NOTE: pressure Pout is assumed to be the uncorrected pressure
!   input from both the OSCAR log and from IPP-generated OPUS headers.
!   pout_corr is the correction defined in the wg_sunrun.dat file

    implicit none
    real(8)::       Pout, xcorr, pout_corr
    character(len=*)::  specname  
    integer(4)::    year, month, day, jday, julian

    read(specname(3:10),'(I4,2I2)')year,month,day
    jday= julian(year,month,day)

!    print *, year,month,day,jday


    if(jday < julian(2008,7,2))then
!       01 July 2008
!       ============
!       Original PTB100 sensor calibrated 01 July 2008
!       See log notes for details, correction = -0.8mb
        Pout = Pout - 0.8

    elseif(jday>julian(2008,7,2) .and. jday<julian(2010,3,4))then

!       New PTB330 sensor installed Dec 2009 and run side-by-side with original PTB100
!       Dec 02 Ronald chnaged calibration by -1 hPa at 1060 hPa
!       Dec 10 revert to old calibration

        if(jday > julian(2009,12,2) .and. jday < julian(2009,12,10)) Pout = Pout - 0.8

!       04 March 2010
!       =============
!       install Vaisala dual sensor PTB330 as main sensor
!       checked against sun photometer sensor, 
!       log records PTB330 from 04 Mar 2010
!       assumed accurate from 4 Mar 2010, no correction necessary (< 0.1 hPa)
!       The old PTB100 reads 0.5-0.65 hPa high at Dec 2009 - Mar 2010
!       Assume linear drift from 0 hPa at 02 July 2008 to +0.6 hPa at 04 Mar 2010

        xcorr = float(jday-julian(2008,7,2))/float(julian(2010,3,4)-julian(2008,7,2))
        Pout = Pout-0.6*xcorr
    else
!       otherwise apply standard pressure correction defined by wg_sunrun.dat
      write(*,*) pout_corr
      pout = pout + pout_corr

    endif

    end subroutine 
