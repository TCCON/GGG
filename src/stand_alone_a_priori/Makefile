# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=stand_alone_a_priori
TSRC  = $(GGGPATH)/src/$(TARGET)
COMN  = $(GGGPATH)/src/comn
GSETUP= $(GGGPATH)/src/gsetup

SOURCES := \
  $(TSRC)/$(TARGET).f \
  $(GSETUP)/calc_aoa.f \
  $(GSETUP)/compute_seasonal_cycle.f \
  $(GSETUP)/apply_secular_trends.f \
  $(GSETUP)/apply_vmr_latitude_gradients.f \
  $(GSETUP)/resample_vmrs_at_effective_altitudes.f \
  $(GSETUP)/read_refvmrs.f \
  $(COMN)/clistindex.f \
  $(COMN)/fbc.f \
  $(COMN)/fnbc.f \
  $(COMN)/lnbc.f \
  $(COMN)/julian.f \
  $(COMN)/lowercase.f \
  $(COMN)/substr.f

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC) $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)
