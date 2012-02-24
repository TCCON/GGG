# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_gfortran

TARGET=create_sunrun_from_Wgong
TARGETDIR=create_sunrun/Wgong

TSRC = $(GGGPATH)/src/$(TARGETDIR)
COMN = $(GGGPATH)/src/comn

SOURCES := \
 $(TSRC)/$(TARGET).f90 \
 $(TSRC)/read_DA8_press.f \
 $(TSRC)/Read_Oscar_log.f90 \
 $(COMN)/read_opus_header.f \
 $(COMN)/getopusparval.f \
 $(COMN)/getendian.f \
 $(COMN)/gindfile.f \
 $(COMN)/lloc.f \
 $(COMN)/lnbc.f \
 $(COMN)/fnbc.f \
 $(COMN)/fbc.f \
 $(COMN)/rbyte.f \
 $(COMN)/write_sunrun.f

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

# Link the object files to create the executable
 $(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Add dependencies to object files
$(TSRC)/create_sunrun_from_Wgong.o: $(GGGPATH)/src/ggg_const_params.f
$(COMN)/gindfile.o:                 $(GGGPATH)/src/ggg_int_params.f
$(COMN)/write_sunrun.o:             $(GGGPATH)/src/ggg_int_params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

