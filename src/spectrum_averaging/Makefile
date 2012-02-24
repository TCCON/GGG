# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=spectrum_averaging
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
  $(TSRC)/$(TARGET).f \
  $(TSRC)/rravgcom.f \
  $(COMN)/read_runlog.f \
  $(COMN)/gindfile.f \
  $(COMN)/lnbc.f \
  $(COMN)/fnbc.f \
  $(COMN)/fbc.f \
  $(COMN)/lloc.f \
  $(COMN)/rbyte.f \
  $(COMN)/substr.f \
  $(COMN)/skiprec.f \
  $(COMN)/getendian.f \
  $(COMN)/write_runlog.f\
  $(COMN)/wlimit.f \
  $(COMN)/get_ggg_environment.f

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Add dependencies to object files
$(TSRC)/rravgcom.o:  $(GGGPATH)/src/ggg_const_params.f $(GGGPATH)/src/ggg_int_params.f
$(COMN)/gindfile.o:  $(GGGPATH)/src/ggg_int_params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

