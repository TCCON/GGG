# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=gsetup
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
 $(TSRC)/$(TARGET).f \
 $(TSRC)/clistindex.f \
 $(TSRC)/write_mav.f \
 $(TSRC)/compute_seasonal_cycle.f \
 $(TSRC)/apply_vmr_latitude_gradients.f \
 $(TSRC)/apply_secular_trends.f \
 $(TSRC)/apply_vmr_latitude_gradients.f \
 $(TSRC)/resample_vmrs_at_effective_altitudes.f \
 $(TSRC)/calc_aoa.f \
 $(TSRC)/readmenu.f \
 $(TSRC)/read_model_fc.f \
 $(TSRC)/read_refvmrs.f \
 $(TSRC)/write_postprocessfile.f \
 $(COMN)/get_ggg_environment.f \
 $(COMN)/fbc.f \
 $(COMN)/fnbc.f \
 $(COMN)/gravity.f \
 $(COMN)/height.f \
 $(COMN)/lnbc.f \
 $(COMN)/vsubs.f \
 $(COMN)/julian.f \
 $(COMN)/caldat.f \
 $(COMN)/lowercase.f \
 $(COMN)/read_runlog.f \
 $(COMN)/read_isotop.f \
 $(COMN)/log1pxox.f \
 $(COMN)/tlpath.f \
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
$(COMN)/gravity.o:                $(GGGPATH)/src/ggg_const_params.f
$(TSRC)/gsetup.o:                 $(GGGPATH)/src/ggg_const_params.f \
                                  $(GGGPATH)/src/ggg_int_params.f $(TSRC)/int_params.f
$(TSRC)/read_model_fc.o:          $(GGGPATH)/src/ggg_const_params.f $(TSRC)/const_params.f
$(TSRC)/write_mav.o:              $(GGGPATH)/src/ggg_const_params.f \
                                  $(GGGPATH)/src/ggg_int_params.f 
$(TSRC)/write_postprocessfile.o:  $(GGGPATH)/src/ggg_int_params.f 
$(COMN)/tlpath.o:                 $(GGGPATH)/src/ggg_const_params.f \
                                  $(GGGPATH)/src/ggg_int_params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

