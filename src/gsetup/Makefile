# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=gsetup
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
 $(TSRC)/$(TARGET).f \
 $(TSRC)/clistindex.f \
 $(TSRC)/write_mav.f \
 $(TSRC)/setvmr.f \
 $(TSRC)/simulate_co2_vmr.f \
 $(TSRC)/gravity.f \
 $(TSRC)/readmenu.f \
 $(TSRC)/readmodFC.f \
 $(TSRC)/readvmrFC.f \
 $(COMN)/fbc.f \
 $(COMN)/fnbc.f \
 $(COMN)/height.f \
 $(COMN)/lnbc.f \
 $(COMN)/lowercase.f \
 $(COMN)/read_runlog.f \
 $(COMN)/read_isotop.f \
 $(COMN)/log1pxox.f \
 $(COMN)/tlpath.f \
 $(COMN)/substr.f 

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#          #  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^
