# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=texo
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
  $(TSRC)/$(TARGET).f \
  $(TSRC)/resample.f \
  $(TSRC)/polyinterp.f \
  $(TSRC)/write_binary_file.f \
  $(COMN)/getendian.f \
  $(COMN)/gindfile.f \
  $(COMN)/gasdev.f \
  $(COMN)/ran1.f \
  $(COMN)/fnbc.f \
  $(COMN)/lloc.f \
  $(COMN)/lnbc.f \
  $(COMN)/profzl.f \
  $(COMN)/rbyte.f \
  $(COMN)/read_runlog.f \
  $(COMN)/write_runlog.f \
  $(COMN)/wlimit.f \
  $(COMN)/vsubs.f
OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

