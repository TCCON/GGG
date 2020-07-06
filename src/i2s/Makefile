include $(GGGPATH)/install/.compiler_ggg

TARGET=i2s

TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
  $(TSRC)/$(TARGET).f \
  $(TSRC)/i2s_processing.f \
  $(TSRC)/get_run_info.f \
  $(TSRC)/check_resistor_info.f \
  $(TSRC)/default_run_info.f \
  $(TSRC)/get_igram_run_parameters.f \
  $(TSRC)/get_opus_xx.f \
  $(TSRC)/get_opusigram_params.f \
  $(TSRC)/get_igram_data.f \
  $(TSRC)/build_cit_name.f \
  $(TSRC)/save_to_file.f \
  $(TSRC)/save_opus.f \
  $(TSRC)/parse_input_top.f \
  $(TSRC)/read_input_line.f \
  $(TSRC)/siv_correction.f \
  $(TSRC)/shars.f \
  $(TSRC)/phase_corr_oper.f \
  $(TSRC)/convolve.f \
  $(TSRC)/resample_ifg.f \
  $(TSRC)/smooth_igram.f \
  $(TSRC)/coeish.f \
  $(COMN)/vrot.f \
  $(COMN)/compute_snr.f \
  $(COMN)/symmetry.f \
  $(COMN)/modfft.f \
  $(COMN)/substr.f \
  $(COMN)/vsubs.f \
  $(COMN)/fnbc.f \
  $(COMN)/fbc.f \
  $(COMN)/lnbc.f \
  $(COMN)/getendian.f \
  $(COMN)/julian.f \
  $(COMN)/rbyte.f \
  $(COMN)/vrot.f \
  $(COMN)/pythag.f \
  $(COMN)/caldat.f \

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

clean:
	rm -f $(OBJECTS) $(TARGET)
