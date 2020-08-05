# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=average_results
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
  $(TSRC)/$(TARGET).f \
  $(TSRC)/covariance_matrix.f \
  $(TSRC)/average_with_mul_bias.f \
  $(TSRC)/average_with_add_bias.f \
  $(COMN)/get_ggg_environment.f \
  $(COMN)/vsubs.f \
  $(COMN)/fbc.f \
  $(COMN)/fnbc.f \
  $(COMN)/lnbc.f \
  $(COMN)/substr.f \
  $(COMN)/read_postproc_header.f \
  $(COMN)/write_postproc_header.f \
  $(COMN)/isclose_s.f \
  $(COMN)/wlimit.f

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Add dependencies to object files
$(TSRC)/average_results.o:        $(GGGPATH)/src/gfit/ggg_int_params.f $(COMN)/postproc_params.f
$(TSRC)/average_with_mul_bias.o:  $(TSRC)/params.f
$(TSRC)/average_with_add_bias.o:  $(TSRC)/params.f
$(COMN)/read_postproc_header.o:   $(COMN)/postproc_params.f
$(COMN)/write_postproc_header.o:  $(COMN)/postproc_params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

