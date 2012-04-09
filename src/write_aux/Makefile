# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=write_aux
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
  $(TSRC)/$(TARGET).f \
  $(TSRC)/read_mav_head.f \
  $(TSRC)/read_mav_aux.f \
  $(COMN)/gravity.f \
  $(COMN)/lnbc.f \
  $(COMN)/skiprec.f

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Add dependencies to object files
$(TSRC)/write_aux.o:  $(GGGPATH)/src/ggg_int_params.f
$(COMN)/gravity.o:  $(GGGPATH)/src/ggg_const_params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

