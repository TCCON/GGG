# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=gsetup
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
 $(TSRC)/$(TARGET).f \
 $(TSRC)/write_mav.f \
 $(TSRC)/read_model_fc.f \
 $(TSRC)/calc_itcz.f \
 $(TSRC)/write_postprocessfile.f \
 $(TSRC)/read_csi.f \
 $(COMN)/zenaz.f \
 $(COMN)/subsolar.f \
 $(COMN)/apply_seasonal_cycle.f \
 $(COMN)/apply_secular_trends.f \
 $(COMN)/apply_vmr_latitude_gradients.f \
 $(COMN)/resample_vmrs_at_effective_altitudes.f \
 $(COMN)/calc_aoa.f \
 $(COMN)/readmenu.f \
 $(COMN)/compute_ztrop.f \
 $(COMN)/compute_ztrop_gct.f \
 $(COMN)/read_refvmrs.f \
 $(COMN)/dmf2wmf.f \
 $(COMN)/clistindex.f \
 $(COMN)/get_ggg_environment.f \
 $(COMN)/fbc.f \
 $(COMN)/fnbc.f \
 $(COMN)/lnbc.f \
 $(COMN)/gravity.f \
 $(COMN)/height.f \
 $(COMN)/lnbc.f \
 $(COMN)/vsubs.f \
 $(COMN)/julian.f \
 $(COMN)/caldat.f \
 $(COMN)/lowercase.f \
 $(COMN)/read_runlog_header.f \
 $(COMN)/read_runlog_data_record.f \
 $(COMN)/read_isotopolog.f \
 $(COMN)/log1pxox.f \
 $(COMN)/tlpath.f \
 $(COMN)/nlo_bs.f \
 $(COMN)/wlimit.f \
 $(COMN)/substr.f 

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Add dependencies to object files
$(TSRC)/gsetup.o:                 $(GGGPATH)/src/gfit/ggg_int_params.f
$(TSRC)/read_model_fc.o:          $(TSRC)/const_params.f
$(TSRC)/write_mav.o:              $(GGGPATH)/src/gfit/ggg_int_params.f 
$(TSRC)/write_postprocessfile.o:  $(GGGPATH)/src/gfit/ggg_int_params.f 
$(COMN)/tlpath.o:                 $(GGGPATH)/src/gfit/ggg_int_params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

