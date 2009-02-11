# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg
TARGET=derive_insitu_correction
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn
CFLAGS = -c -O -fbounds-check -Wunused -Wuninitialized -W -Wall
LDFLAGS =

#  Statement of dependencies
$(TARGET) : $(TSRC)/$(TARGET).o \
  $(COMN)/lnbc.f \
  $(COMN)/fnbc.f \
  $(COMN)/fbc.f \
  $(COMN)/substr.f

#  Link the object files to create the executable
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Compile the .f files that are newer than the corresponding .o files
%.o : %.f
	$(FC)  $(CFLAGS) -o $@ $<
