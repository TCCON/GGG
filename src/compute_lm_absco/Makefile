# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg
ifeq ($(FC),gfortran)
      LOCAL_CFLAGS=-ffloat-store
else
      LOCAL_CFLAGS=
endif

TARGET=compute_lm_absco
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
 $(TSRC)/$(TARGET).f \
 $(TSRC)/o2_lm_cia_routines.f \
 $(TSRC)/co2_all_lm_routines.f \
 $(TSRC)/ch4_all_lm_routines.f \
 $(TSRC)/humlik.f \
 $(TSRC)/cpf.f \
 $(COMN)/get_ggg_environment.f \
 $(COMN)/fbc.f \
 $(COMN)/fnbc.f \
 $(COMN)/lnbc.f \
 $(COMN)/substr.f 

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) $(LOCAL_CFLAGS) -o $@ $<

#          #  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Add dependencies to object files
$(TSRC)/compute_lm_absco.o:  $(GGGPATH)/src/gfit/ggg_int_params.f
$(TSRC)/o2_lm_cia_routines.o:  $(GGGPATH)/src/gfit/ggg_int_params.f
$(TSRC)/co2_all_lm_routines.o:  $(GGGPATH)/src/gfit/ggg_int_params.f
$(TSRC)/ch4_all_lm_routines.o:  $(GGGPATH)/src/gfit/ggg_int_params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)
