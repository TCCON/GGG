# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=bin2asc
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
  $(TSRC)/$(TARGET).f \
  $(COMN)/gindfile.f \
  $(COMN)/lnbc.f \
  $(COMN)/lloc.f \
  $(COMN)/rbyte.f \
  $(COMN)/getendian.f \
  $(COMN)/write_runlog.f \
  $(COMN)/wlimit.f \
  $(COMN)/read_runlog.f \
  $(COMN)/get_ggg_environment.f

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Add dependencies to object files
$(TSRC)/bin2asc.o:  $(GGGPATH)/src/ggg_int_params.f
$(COMN)/gindfile.o: $(GGGPATH)/src/ggg_int_params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

