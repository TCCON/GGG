#include $(GGGPATH)/install/.compiler_ggg
FC=gfortran
#FC=ifort
CFLAGS = -c -O -fbounds-check -Wuninitialized -Wall -W -g -static
#CFLAGS = -c -assume byterecl -fpe0 -g -traceback -check bounds -heap-arrays 32
LDFLAGS =

TARGET=opus-i2s

SOURCES := \
  $(TARGET).f \
  ../comn/parse_input_top.f \
  ../comn/i2s_lite_siv.f \
  ../comn/siv_correction.f \
  ../comn/shars.f \
  ../comn/resample_ifg.f \
  ../comn/smooth_igram.f \
  ../comn/coeish.f \
  ../comn/vrot.f \
  ../comn/modfft.f \
  ../comn/substr.f \
  ../comn/vsubs.f \
  ../comn/fnbc.f \
  ../comn/fbc.f \
  ../comn/lnbc.f \
  ../comn/getendian.f \
  ../comn/julian.f \
  ../comn/xxrev.f \
  ../comn/coeish.f \
  ../comn/vrot.f \
  ../comn/pythag.f \
  ../comn/resample_ifg.f \
  ../opus-comn/read_input_line.f \
  ../opus-comn/get_opus_xx.f \
  ../opus-comn/get_opusigram_params.f \
  ../opus-comn/build_cit_name.f \
  ../opus-comn/save_to_file.f

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $@ $^

clean:
	rm -f $(OBJECTS) $(TARGET)

