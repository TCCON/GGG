# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

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
 $(COMN)/get_ggg_environment.f \
 $(COMN)/write_sunrun.f

OBJECTSF := $(patsubst %.f, %.o, $(filter %.f, $(SOURCES)))
OBJECTSF90 := $(patsubst %.f90, %.o, $(filter %.f90, $(SOURCES)))
OBJECTS := $(OBJECTSF) $(OBJECTSF90)

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Compile the .f90 files that are newer than the corresponding .o files
%.o: %.f90
	$(FC)  $(CFLAGS) -o $@ $<

# Link the object files to create the executable
 $(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Add dependencies to object files
$(COMN)/gindfile.o:                 $(GGGPATH)/src/gfit/ggg_int_params.f
$(COMN)/write_sunrun.o:             $(GGGPATH)/src/gfit/ggg_int_params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

