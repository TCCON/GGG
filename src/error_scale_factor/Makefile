# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=error_scale_factor
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
  $(TSRC)/$(TARGET).f \
  $(COMN)/lnbc.f \
  $(COMN)/lowercase.f \
  $(COMN)/fnbc.f \
  $(COMN)/fbc.f \
  $(COMN)/substr.f \
  $(COMN)/read_postproc_header.f

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Add necessary parameter prerequisite files
$(TSRC)/$(TARGET).o:             $(COMN)/postproc_params.f
$(COMN)/read_postproc_header.o:  $(COMN)/postproc_params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

