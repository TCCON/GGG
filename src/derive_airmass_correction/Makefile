# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=derive_airmass_correction
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
  $(TSRC)/$(TARGET).f \
  $(COMN)/wlsfit.f \
  $(COMN)/hfti.f \
  $(COMN)/lnbc.f \
  $(COMN)/fnbc.f \
  $(COMN)/fbc.f \
  $(COMN)/lloc.f \
  $(COMN)/lowercase.f \
  $(COMN)/substr.f

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Add dependencies to object files
$(TSRC)/derive_airmass_correction.o: $(GGGPATH)/src/gfit/ggg_int_params.f $(COMN)/postproc_params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

