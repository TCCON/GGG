# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=write_eof
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
  $(TSRC)/$(TARGET).f \
  $(TSRC)/prepare_oof_output.f \
  $(TSRC)/prepare_collate_all.f \
  $(TSRC)/read_oneline_oof.f \
  $(TSRC)/read_oneline_col.f \
  $(COMN)/lnbc.f \
  $(COMN)/fndlun.f \
  $(COMN)/getlun.f \
  $(COMN)/fnbc.f \
  $(COMN)/read_runlog.f \
  $(COMN)/lloc.f \
  $(COMN)/substr.f \
  $(COMN)/fnbc.f \
  $(COMN)/fbc.f \
  $(COMN)/wlimit.f \
  $(COMN)/skiprec.f

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^


clean:
	rm -f *.o