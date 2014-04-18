# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=generate_spts
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
 $(TSRC)/$(TARGET).f \
 $(COMN)/lnbc.f \
 $(COMN)/file_size_in_bytes.f \
 $(COMN)/solar_pts.f \
 $(COMN)/get_ggg_environment.f \
 $(COMN)/posnall.f 

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

