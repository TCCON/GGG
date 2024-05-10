# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=apply_tcor
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
 $(TSRC)/$(TARGET).f \
 $(TSRC)/compute_tdel.f \
 $(COMN)/clistindex.f \
 $(COMN)/get_ggg_environment.f \
 $(COMN)/fbc.f \
 $(COMN)/fnbc.f \
 $(COMN)/lnbc.f \
 $(COMN)/vsubs.f \
 $(COMN)/substr.f 

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Add dependencies to object files
$(TSRC)/$(TARGET).o:       $(COMN)/postproc_params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

