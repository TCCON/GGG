FC = g77
TARGET=create_runlog
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn
CFLAGS = -c -O -fbounds-check -Wunused -Wuninitialized -W -Wall
LDFLAGS =

#  Statement of dependencies
$(TARGET): \
 $(TSRC)/$(TARGET).o \
 $(TSRC)/rdsphead.o \
 $(COMN)/read_opus_header.o \
 $(COMN)/getopusparval.o \
 $(COMN)/read_mkiv_header.o \
 $(COMN)/getendian.o \
 $(COMN)/gindfile.o \
 $(COMN)/lloc.o \
 $(COMN)/julian.o \
 $(COMN)/lnbc.o \
 $(COMN)/fnbc.o \
 $(COMN)/fbc.o \
 $(COMN)/rbyte.o \
 $(COMN)/subsolar.o \
 $(COMN)/read_sunrun.o \
 $(COMN)/write_runlog.o \
 $(COMN)/wlimit.o \
 $(COMN)/zenaz.o

#  Link the object files to create the executable
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Compile the .f files that are newer than the corresponding .o files
%.o : %.f
	$(FC)  $(CFLAGS) -o $@ $<
