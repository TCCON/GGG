# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET= create_sunrun_from_bremen

TSRC = $(GGGPATH)/src/create_sunrun/bremen
COMN = $(GGGPATH)/src/comn

SOURCES := \
 $(TSRC)/$(TARGET).f \
 $(TSRC)/get_bremen_log_vals.f \
 $(COMN)/read_opus_header.f \
 $(COMN)/getopusparval.f \
 $(COMN)/getendian.f \
 $(COMN)/get_ggg_environment.f \
 $(COMN)/gindfile.f \
 $(COMN)/lloc.f \
 $(COMN)/lnbc.f \
 $(COMN)/fnbc.f \
 $(COMN)/fbc.f \
 $(COMN)/julian.f \
 $(COMN)/substr.f \
 $(COMN)/rbyte.f \
 $(COMN)/write_sunrun.f 

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#          #  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

