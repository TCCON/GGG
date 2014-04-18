# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=diurnret
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
  $(TSRC)/$(TARGET).f \
  $(TSRC)/newnnls.f \
  $(TSRC)/nnls.f \
  $(COMN)/get_ggg_environment.f \
  $(COMN)/lowercase.f \
  $(COMN)/subsolar.f \
  $(COMN)/get_ggg_environment.f \
  $(COMN)/twodint.f \
  $(COMN)/zenaz.f \
  $(COMN)/julian.f \
  $(COMN)/fbc.f \
  $(COMN)/fnbc.f \
  $(COMN)/lnbc.f \
  $(COMN)/substr.f

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Add dependencies to object files
$(TSRC)/diurnret.o:  $(GGGPATH)/src/ggg_int_params.f $(GGGPATH)/src/ggg_const_params.f
$(COMN)/zenaz.o:     $(GGGPATH)/src/ggg_const_params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

