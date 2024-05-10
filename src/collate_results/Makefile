# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=collate_results
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
  $(TSRC)/$(TARGET).f \
  $(TSRC)/generate_qc_flag.f \
  $(COMN)/fbc.f \
  $(COMN)/fnbc.f \
  $(COMN)/lnbc.f \
  $(COMN)/lloc.f \
  $(COMN)/read_runlog_header.f \
  $(COMN)/read_runlog_data_record.f \
  $(COMN)/get_ggg_environment.f \
  $(COMN)/substr.f \
  $(COMN)/write_postproc_header.f \
  $(COMN)/isclose_s.f

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Add dependencies to object files
$(TSRC)/collate_results.o:  $(GGGPATH)/src/gfit/ggg_int_params.f $(TSRC)/params.f $(COMN)/postproc_params.f
$(TSRC)/generate_qc_flag.o: $(TSRC)/params.f
$(COMN)/write_postproc_header.o: $(COMN)/postproc_params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

