# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=convert_runlog
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
 $(TSRC)/$(TARGET).f \
 $(COMN)/read_runlog_data_record.f \
 $(COMN)/read_runlog_header.f \
 $(COMN)/write_runlog_header.f \
 $(COMN)/write_runlog_data_record.f \
 $(COMN)/wlimit.f \
 $(COMN)/substr.f \
 $(COMN)/lnbc.f \
 $(COMN)/fnbc.f \
 $(COMN)/fbc.f

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#          #  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

