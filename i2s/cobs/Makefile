#include $(GGGPATH)/install/.compiler_ggg
FC=gfortran
CFLAGS = -c -O -fbounds-check -Wuninitialized -Wall -W
#FC=ifort
#CFLAGS = -c -assume byterecl -fpe0 -g -traceback -check bounds -heap-arrays 32
LDFLAGS =

TARGET=cobs

SOURCES := \
  $(TARGET).f 

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $@ $^

clean:
	rm -f $(OBJECTS) $(TARGET)

