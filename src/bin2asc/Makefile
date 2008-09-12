FC = g77
TARGET=bin2asc
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn
CFLAGS = -c -O -fbounds-check -Wunused -Wuninitialized -W -Wall
LDFLAGS =

#  Statement of dependencies
$(TARGET) : $(TSRC)/$(TARGET).o \
  $(COMN)/gindfile.o \
  $(COMN)/lnbc.o \
  $(COMN)/lloc.o \
  $(COMN)/rbyte.o \
  $(COMN)/getendian.o \
  $(COMN)/write_runlog.o\
  $(COMN)/wlimit.o\
  $(COMN)/read_runlog.o

#  Link the object files to create the executable
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Compile the .f files that are newer than the corresponding .o files
%.o : %.f
	$(FC)  $(CFLAGS) -o $@ $<
