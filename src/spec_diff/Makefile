# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=spec_diff
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
 $(TSRC)/$(TARGET).f \
 $(COMN)/rdsphead.f \
 $(COMN)/read_DA8_header.f \
 $(COMN)/read_opus_header.f \
 $(COMN)/getopusparval.f \
 $(COMN)/read_mkiv_header.f \
 $(COMN)/getendian.f \
 $(COMN)/gindfile.f \
 $(COMN)/lloc.f \
 $(COMN)/julian.f \
 $(COMN)/lnbc.f \
 $(COMN)/fnbc.f \
 $(COMN)/fbc.f \
 $(COMN)/rbyte.f \
 $(COMN)/subsolar.f \
 $(COMN)/substr.f \
 $(COMN)/read_sunrun.f \
 $(COMN)/write_runlog_header.f \
 $(COMN)/write_runlog_data_record.f \
 $(COMN)/wlimit.f \
 $(COMN)/zenaz.f \
 $(COMN)/get_ggg_environment.f

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Add dependencies to object files
$(TSRC)/create_runlog.o:  $(GGGPATH)/src/gfit/ggg_int_params.f
$(COMN)/gindfile.o:       $(GGGPATH)/src/gfit/ggg_int_params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

