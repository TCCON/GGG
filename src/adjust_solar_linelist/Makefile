# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=adjust_solar_linelist
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn
#CFLAGS = -c -O -fbounds-check -Wunused -Wuninitialized -W -Wall
#LDFLAGS =

#  Statement of dependencies
SOURCES := \
        $(TSRC)/$(TARGET).f \
	$(COMN)/fbc.f \
	$(COMN)/fnbc.f \
	$(COMN)/lnbc.f \
	$(COMN)/posnall.f \
	$(COMN)/file_size_in_bytes.f \
	$(COMN)/vsubs.f \
	$(COMN)/hfti.f

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -mcmodel=medium -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Add dependencies to object files
$(TSRC)/adjust_solar_linelist.o:  $(TSRC)/params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

