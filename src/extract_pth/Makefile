# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=extract_pth
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
  $(TSRC)/$(TARGET).f \
  $(COMN)/log1pxox.o \
  $(COMN)/get_ggg_environment.o \
  $(COMN)/read_runlog_header.o \
  $(COMN)/read_runlog_data_record.o \
  $(COMN)/substr.o \
  $(COMN)/julian.o \
  $(COMN)/caldat.o \
  $(COMN)/fbc.o \
  $(COMN)/fnbc.o \
  $(COMN)/lnbc.o

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

