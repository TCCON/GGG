# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET= create_sunrun

TSRC = $(GGGPATH)/src/create_sunrun
COMN = $(GGGPATH)/src/comn

SOURCES := \
 $(TSRC)/$(TARGET).f \
 $(COMN)/read_opus_header.f \
 $(COMN)/read_mkiv_header.f \
 $(COMN)/get_ggg_environment.f \
 $(COMN)/getopusparval.f \
 $(COMN)/getendian.f \
 $(COMN)/gindfile.f \
 $(COMN)/lloc.f \
 $(COMN)/lnbc.f \
 $(COMN)/fnbc.f \
 $(COMN)/fbc.f \
 $(COMN)/rbyte.f \
 $(COMN)/substr.f \
 $(COMN)/write_sunrun.f \
 $(COMN)/get_ggg_environment.f

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Add dependencies to object files
$(TSRC)/create_sunrun.o:  $(GGGPATH)/src/ggg_int_params.f $(TSRC)/params.f
$(COMN)/gindfile.o:       $(GGGPATH)/src/ggg_int_params.f
$(COMN)/write_sunrun.o:   $(GGGPATH)/src/ggg_int_params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

