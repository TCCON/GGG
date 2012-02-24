# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=apply_airmass_correction
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
  $(TSRC)/$(TARGET).f \
  $(COMN)/lnbc.f \
  $(COMN)/fnbc.f \
  $(COMN)/fbc.f \
  $(COMN)/substr.f \
  $(COMN)/get_ggg_environment.f

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Add dependencies to object files
$(TSRC)/apply_airmass_correction.o:  $(GGGPATH)/src/ggg_int_params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

