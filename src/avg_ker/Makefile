FC = g77
TARGET=avg_ker
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn
CFLAGS = -c -O -fbounds-check -Wunused -Wuninitialized -W -Wall
LDFLAGS =

#  Statement of dependencies
$(TARGET) : $(TSRC)/$(TARGET).o \
	$(COMN)/fbc.f \
	$(COMN)/fnbc.f \
	$(COMN)/lnbc.f \
	$(COMN)/vsubs.f \
	$(COMN)/hfti.f

#  Link the object files to create the executable
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Compile the .f files that are newer than the corresponding .o files
%.o : %.f
	$(FC)  $(CFLAGS) -o $@ $<
