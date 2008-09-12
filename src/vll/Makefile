FC = g77
TARGET=vll
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn
CFLAGS = -c -O -fbounds-check -Wunused -Wuninitialized -W -Wall
LDFLAGS =

#  Statement of dependencies
$(TARGET) : $(TSRC)/$(TARGET).o \
  $(COMN)/lnbc.o \
  $(COMN)/lloc.o \
  $(COMN)/hitran_to_atmos_gas_numbering.o \
  $(COMN)/file_size_in_bytes.o \
  $(COMN)/posnall.o

#  Link the object files to create the executable
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Compile the .f files that are newer than the corresponding .o files
%.o : %.f
	$(FC)  $(CFLAGS) -o $@ $<
