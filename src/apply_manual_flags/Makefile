# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=apply_manual_flags
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
  $(TSRC)/$(TARGET).f \
  $(COMN)/lnbc.f \
  $(COMN)/fnbc.f \
  $(COMN)/fbc.f \
  $(COMN)/substr.f \
  $(COMN)/get_ggg_environment.f \
  $(COMN)/read_postproc_header.f \
  $(COMN)/write_postproc_header.f \
  $(COMN)/isclose_d.f \
  $(COMN)/julian.f

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Add dependencies to object files
$(TSRC)/apply_insitu_correction.o:  $(GGGPATH)/src/gfit/ggg_int_params.f $(COMN)/postproc_params.f
$(COMN)/read_postproc_header.o:     $(COMN)/postproc_params.f
$(COMN)/write_postproc_header.o:     $(COMN)/postproc_params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)