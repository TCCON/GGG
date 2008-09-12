FC = g77
TARGET=gggavg
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn
CFLAGS = -c -O -fbounds-check -Wunused -Wuninitialized -W -Wall
LDFLAGS =

#  Statement of dependencies
$(TARGET): \
  $(TSRC)/$(TARGET).o \
  $(TSRC)/t_fit.o \
  $(COMN)/fbc.o \
  $(COMN)/fnbc.o \
  $(COMN)/lnbc.o \
  $(COMN)/read_runlog.o \
  $(COMN)/substr.o

#  Link the object files to create the executable
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Compile the .f files that are newer than the corresponding .o files
%.o : %.f
	$(FC)  $(CFLAGS) -o $@ $<

