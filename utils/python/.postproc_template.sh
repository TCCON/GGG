$GGGPATH/bin/collate_results t
$GGGPATH/bin/collate_results v
$GGGPATH/bin/average_results $RUNLOG.tsw
$GGGPATH/bin/average_results $RUNLOG.vsw
$GGGPATH/bin/apply_airmass_correction $RUNLOG.vsw
$GGGPATH/bin/average_results $RUNLOG.vsw.ada
$GGGPATH/bin/apply_insitu_correction $RUNLOG.vav.ada
$GGGPATH/bin/error_scale_factor $RUNLOG.vav.ada.aia
$GGGPATH/bin/extract_pth $RUNLOG.grl y
$GGGPATH/bin/write_official_output_file $RUNLOG.vav.ada.aia
$GGGPATH/bin/apply_manual_flags $RUNLOG.vav.ada.aia.oof
$GGGPATH/bin/write_netcdf $RUNLOG.tav
$GGGPATH/bin/write_aux $RUNLOG.mav
